#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Run the multiple alignment on pairwise alignment input (ie cactus_setup_phase and beyond)
   
"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.cactus_progressive import setupBinaries, importSingularityImage
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.pipeline.cactus_workflow import CactusSetupCheckpoint
from cactus.shared.common import makeURL


from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("blastOutput", type=str, help = "Blast output (from cactus-blast)")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--database", dest="database",
                      help="Database type: tokyo_cabinet or kyoto_tycoon"
                      " [default: %(default)s]",
                      default="kyoto_tycoon")
    parser.add_argument("--configFile", dest="configFile",
                      help="Specify cactus configuration file",
                      default=None)
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                        " must appear in NEWICK tree in <seqfile>) to use as a "
                        "root for the alignment.  Any genomes not below this node "
                        "in the tree may be used as outgroups but will never appear"
                        " in the output.  If no root is specifed then the root"
                        " of the tree is used. ", default=None, required=True)
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    setLoggingFromOptions(options)

    options.database = 'kyoto_tycoon'
    
    # Mess with some toil options to create useful defaults.

    # Caching generally slows down the cactus workflow, plus some
    # methods like readGlobalFileStream don't support forced
    # reads directly from the job store rather than from cache.
    options.disableCaching = True
    # Job chaining breaks service termination timing, causing unused
    # databases to accumulate and waste memory for no reason.
    options.disableChaining = True
    if options.retryCount is None:
        # If the user didn't specify a retryCount value, make it 5
        # instead of Toil's default (1).
        options.retryCount = 5

    start_time = timeit.default_timer()
    runCactusBlastOnly(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-blast has finished after {} seconds".format(run_time))

def runCactusBlastOnly(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            alignmentID = toil.restart()
        else:
            options.cactusDir = getTempDirectory()
            
            #Create the progressive cactus project (as we do in runCactusProgressive)
            projWrapper = ProjectWrapper(options)
            projWrapper.writeXml()

            pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                                  '%s_project.xml' % ProjectWrapper.alignmentDirName)
            assert os.path.exists(pjPath)

            project = MultiCactusProject()

            if not os.path.isdir(options.cactusDir):
                os.makedirs(options.cactusDir)

            project.readXML(pjPath)

            # open up the experiment (as we do in ProgressiveUp.run)
            # note that we copy the path into the options here
            experimentFile = project.expMap[options.root]
            expXml = ET.parse(experimentFile).getroot()
            experiment = ExperimentWrapper(expXml)
            configPath = experiment.getConfigPath()
            configXml = ET.parse(configPath).getroot()
            
            seqIDMap = dict()
            tree = MultiCactusTree(experiment.getTree()).extractSubTree(options.root)
            leaves = [tree.getName(leaf) for leaf in tree.getLeaves()]
            outgroups = experiment.getOutgroupGenomes()
            genome_set = set(leaves + outgroups)

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            for genome, seq in list(project.inputSequenceMap.items()):
                if genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    experiment.setSequenceID(genome, toil.importFile(seq))

            # import the outgroups
            outgroupIDs = []
            for i, outgroup in enumerate(outgroups):
                outgroupID = toil.importFile(makeURL(options.blastOutput) + '.og_fragment_{}'.format(i))
                outgroupIDs.append(outgroupID)
                experiment.setSequenceID(outgroup, outgroupID)

            # write back the experiment, as CactusWorkflowArguments wants a path
            experiment.writeXML(experimentFile)

            #import cactus config
            if options.configFile:
                cactusConfigID = toil.importFile(makeURL(options.configFile))
            else:
                cactusConfigID = toil.importFile(makeURL(project.getConfigPath()))
            project.setConfigID(cactusConfigID)

            project.syncToFileStore(toil)
            configNode = ET.parse(project.getConfigPath()).getroot()
            configWrapper = ConfigWrapper(configNode)
            configWrapper.substituteAllPredefinedConstantsWithLiterals()            
            
            workFlowArgs = CactusWorkflowArguments(options, experimentFile=experimentFile, configNode=configNode, seqIDMap = project.inputSequenceIDMap)

            # hack
            totalSequenceSizeID = toil.importFile(makeURL(options.blastOutput) + '.total_sequence_size')

            #import the files that cactus-blast made
            workFlowArgs.alignmentsID = toil.importFile(makeURL(options.blastOutput))
            try:
                workFlowArgs.secondaryAlignmentsID = toil.importFile(makeURL(options.blastOutput) + '.secondary')
            except:
                workFlowArgs.secondaryAlignmentsID = None
            workFlowArgs.outgroupFragmentIDs = outgroupIDs
            workFlowArgs.ingroupCoverageIDs = []
            for i in range(len(leaves)):
                workFlowArgs.ingroupCoverageIDs.append(toil.importFile(makeURL(options.blastOutput) + '.ig_coverage_{}'.format(i)))

            halID = toil.start(Job.wrapJobFn(run_setup_phase, workFlowArgs, totalSequenceSizeID))

        # export the hal
        toil.exportFile(halID, makeURL(options.outputHal))

def run_setup_phase(job, cactusWorkflowArguments, totalSequenceSizeID):
    """ hack to get the sequence size from a file into the arguments """
    totalSizePath = os.path.join(job.fileStore.getLocalTempDir(), 'totalSequenceSize')
    job.fileStore.readGlobalFile(totalSequenceSizeID, totalSizePath)
    with open(totalSizePath, 'r') as sizeFile:
        cactusWorkflowArguments.totalSequenceSize = int(sizeFile.read().strip())
    return job.addChild(CactusSetupCheckpoint(cactusWorkflowArguments=cactusWorkflowArguments, phaseName="setup")).rv()

if __name__ == '__main__':
    main()
