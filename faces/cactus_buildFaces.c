#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "cactusGlobalsPrivate.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"

#include "cactus_buildFaces.h"

typedef struct _liftedEdge LiftedEdge;



////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Ancestors function
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Utility function for the lifted edge hashtable
 */
static uint32_t buildFaces_hashfunction(void *ptr) {
	Cap *key = (Cap *) ptr;
	return (uint32_t) cap_getName(key);
}

/*
 * Utility function for the lifted edge hashtable
 */
static int32_t buildFaces_key_eq_fn(void *ptrA, void *ptrB) {
	return ptrA == ptrB;
}

/*
 * compute hash table of ancestors
 */
static struct hashtable * buildFaces_computeAncestors(Net * net) {
	struct hashtable *ancestorsTable = create_hashtable(16,
			buildFaces_hashfunction, buildFaces_key_eq_fn, NULL, NULL);
	Net_CapIterator *iter = net_getCapIterator(net);
	Cap * cap, * tmp;

	logInfo("Computing ancestors\n");

	while ((cap = net_getNextCap(iter))) {
		// ... check if connected
		if (cap_getAdjacency(cap) && cap_getParent(cap)) {
			// Go up the tree
			for (tmp = cap_getParent(cap); tmp; tmp = cap_getParent(tmp)) {
				if (cap_getAdjacency(tmp) || !cap_getParent(tmp)) {
					hashtable_insert(ancestorsTable, cap, tmp);
					break;
				}
			}
		}
	}

	net_destructCapIterator(iter);
	return ancestorsTable;

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Lifted edges functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _liftedEdge {
	Cap *destination;
	Cap *bottomNode;
};

/*
 * Lifted edge destructor
 */
static void buildFaces_destructLiftedEdge(LiftedEdge * liftedEdge) {
	free(liftedEdge);
}
/*
 * Utility function for the lifted edge hashtable
 */
static void buildFaces_destructValue(void *ptr) {
	destructList((struct List *) ptr);
}

/*
 * Utility function for List struct
 */
static void buildFaces_destructListElem(void *ptr) {
	buildFaces_destructLiftedEdge((LiftedEdge *) ptr);
}

/*
 * Fill in a hashtable which to every node associates
 * a list of lifted edges
 */
static struct hashtable *buildFaces_computeLiftedEdges(Net * net, struct hashtable * ancestorsTable) {
	struct hashtable *liftedEdgesTable = create_hashtable(16,
			buildFaces_hashfunction, buildFaces_key_eq_fn, NULL,
			buildFaces_destructValue);
	Net_CapIterator *iter = net_getCapIterator(net);
	Cap *cap, *attachedAncestor;
	Cap *adjacency, *adjacencyAncestor;
	struct List *liftedEdges;
	LiftedEdge *liftedEdge;

	logInfo("Computing lifted edges\n");

	// Iterate through potential bottom nodes
	while ((cap = net_getNextCap(iter))) {
		// ... check if connected
		if ((adjacency = cap_getAdjacency(cap))) {
			// ... lift
			attachedAncestor = hashtable_search(ancestorsTable, cap);
			adjacencyAncestor = hashtable_search(ancestorsTable, cap_getPositiveOrientation(adjacency));

#ifdef BEN_DEBUG
			assert((attachedAncestor && adjacencyAncestor) || (!attachedAncestor && !adjacencyAncestor));
#endif 

			// If root node
			if (attachedAncestor == NULL)
				continue;

			// ... create lifted edge
			liftedEdge = mallocLocal(sizeof(LiftedEdge));
			liftedEdge->destination = adjacencyAncestor;
			liftedEdge->bottomNode = cap;

#ifdef BEN_DEBUG
			// Self loop
			if (adjacencyAncestor == attachedAncestor)
				abort();
#endif

			// ... add it to the hashtable
			if ((liftedEdges = hashtable_search(liftedEdgesTable,
					attachedAncestor))) {
				listAppend(liftedEdges, liftedEdge);
			} else {
				liftedEdges = constructZeroLengthList(2,
						buildFaces_destructListElem);
				listAppend(liftedEdges, liftedEdge);
				hashtable_insert(liftedEdgesTable, attachedAncestor,
						liftedEdges);
			}
		}
	}

	net_destructCapIterator(iter);
	return liftedEdgesTable;
}

/*
 * Recursive function which fills a given list with the
 * connected nodes within a module
 */
static void buildFaces_fillTopNodeList(Cap * cap, struct List *list,
		struct hashtable *liftedEdgesTable) {
	struct List *liftedEdges;
	int32_t index;

	// Limit of recursion
	if (listContains(list, cap))
		return;

	// Actual filling
	logInfo("Adding cap %p to face\n", cap);
	listAppend(list, cap);

	// Recursion through lifted edges
	if ((liftedEdges = hashtable_search(liftedEdgesTable, cap)))
		for (index = 0; index < liftedEdges->length; index++)
			buildFaces_fillTopNodeList(
					((LiftedEdge *) liftedEdges-> list[index])->destination,
					list, liftedEdgesTable);

	// Recursion through adjacency
	if (cap_getAdjacency(cap))
		buildFaces_fillTopNodeList(cap_getAdjacency(cap), list,
				liftedEdgesTable);

	// Remove from lifted edges table to prevent double usage
	// of end instances
	hashtable_remove(liftedEdgesTable, cap, 0);
	if (liftedEdges)
		destructList(liftedEdges);
}

/*
 * Produces the end instance which is the destination of the unique lifted edge
 * out of a top node
 */
static Cap *buildFaces_getMinorLiftedEdgeDestination(Cap * cap,
		struct List *liftedEdges) {
	int32_t index;
	Cap * adjacency = cap_getAdjacency(cap);
	Cap *ancestralEdgeDestination = NULL;
	Cap *liftedDestination;

	if (adjacency)
		ancestralEdgeDestination = cap_getPositiveOrientation(cap_getAdjacency(
				cap));

#ifndef BEN_DEBUG
	for (index = 0; index < liftedEdges->length; index++)
	if ((liftedDestination = ((LiftedEdge *) liftedEdges->list[index])->destination)
			&& liftedDestination != ancestralEdgeDestination)
	return liftedDestination;

	return NULL;
#else
	Cap * candidate = NULL;

	// Ensure that no more than one derived edge destinations
	for (index = 0; index < liftedEdges->length; index++) {
		if ((liftedDestination
				= ((LiftedEdge *) liftedEdges->list[index])->destination)
				&& liftedDestination != ancestralEdgeDestination) {
			if (!candidate)
				candidate = liftedDestination;
			else
				abort();
		}
	}
	return candidate;
#endif
}

/*
 * Constructs a face from a given Cap
 */
static void buildFaces_constructFromCap(Cap * startingCap,
		struct hashtable *liftedEdgesTable, struct hashtable *ancestorsTable, Net * net) {
	Face *face = face_construct(net);
	struct List *topNodes = constructZeroLengthList(16, NULL);
	struct List *liftedEdges;
	Cap *cap, *bottomNode, *ancestor;
	int32_t index, index2;

	printf("Constructing new face");

	// Establish list of top nodes
	buildFaces_fillTopNodeList(startingCap, topNodes, liftedEdgesTable);

#ifdef BEN_DEBUG
	// What, no top nodes!?
	if (topNodes->length == 0)
		abort();
#endif

	printf("Cardinal = %i\n", topNodes->length);

	// Initialize data structure
	face_allocateSpace(face, topNodes->length);

	// For every top node
	for (index = 0; index < topNodes->length; index++) {
		cap = topNodes->list[index];
		face_setTopNode(face, index, cap);
		liftedEdges = hashtable_search(liftedEdgesTable, cap);

		if (!liftedEdges) {
			face_setBottomNodeNumber(face, index, 0);
			continue;
		}

		face_setBottomNodeNumber(face, index, liftedEdges->length);
		// For every bottom node of that top node
		for (index2 = 0; index2 < liftedEdges->length; index2++) {
			bottomNode = ((LiftedEdge *) liftedEdges-> list[index2])->bottomNode;
			face_addBottomNode(face, index, bottomNode);

#if BEN_DEBUG
			assert(cap_getAdjacency(bottomNode));
#endif 
			ancestor = hashtable_search(ancestorsTable, cap_getPositiveOrientation(cap_getAdjacency(bottomNode)));
			if (cap_getAdjacency(cap) != ancestor) 
				face_setDerivedDestination(face, index, index2, ancestor);
			else
				face_setDerivedDestination(face, index, index2, NULL);
				
#ifdef BEN_DEBUG
			// If bottom nodes part of top nodes
			if (listContains(topNodes, cap_getPositiveOrientation(
					((LiftedEdge*) liftedEdges->list[index2])->bottomNode)))
				abort();
#endif
		}
	}

#ifdef BEN_DEBUG_ULTRA
	if (!buildFaces_isSimple(face))
	abort();
#endif

	// Clean up
	destructList(topNodes);
}

/*
 * Construct faces in net and add them to the Net's pointers
 */
void buildFaces_constructFaces(Net * net) {
	struct hashtable * ancestorsTable = buildFaces_computeAncestors(net);
	struct hashtable *liftedEdgesTable = buildFaces_computeLiftedEdges(net, ancestorsTable);
	Net_CapIterator *iter = net_getCapIterator(net);
	struct List *liftedEdges;
	Cap *current;

	logInfo("Constructing faces\n");

	while ((current = net_getNextCap(iter)))
		if ((liftedEdges = hashtable_search(liftedEdgesTable, current))
		    && buildFaces_getMinorLiftedEdgeDestination(current, liftedEdges))
			buildFaces_constructFromCap(current, liftedEdgesTable, ancestorsTable, net);

	hashtable_destroy(liftedEdgesTable, true,false);
	hashtable_destroy(ancestorsTable, false,false);
	net_destructCapIterator(iter);
}

/*
 * Simplify a given face	
 */
void buildFaces_simplify(Face * face, Net * net) {
	// TODO
	// Here lies NP-completeness

#ifdef BEN_DEBUG
	assert(face_isSimple(face));
#endif 

}

/*
 * Simplify all the faces in the net
 */
void buildFaces_simplifyFaces(Net * net) {
	Net_FaceIterator *iter = net_getFaceIterator(net);
	Face *face;

	logInfo("Simplifying faces\n");

	while ((face = net_getNextFace(iter)))
		buildFaces_simplify(face, net);

	net_destructFaceIterator(iter);
}

/*
 * Create an interpolation between parent and child caps at event
 */
static Cap *buildFaces_interpolateCaps(Cap * parentCap, Cap * childCap,
		Event * event) {
	Cap * newCap = cap_construct(cap_getEnd(childCap), event);
	cap_changeParentAndChild(newCap, childCap);
	cap_makeParentAndChild(parentCap, newCap);
	return newCap;
}

/*
 * Creates an inteprolation half way on the branch between two events
 */
static Event *buildFaces_interpolateEvents(Event* parentEvent,
		Event* childEvent) {
	EventTree * eventTree = event_getEventTree(parentEvent);
	float branchLength = event_getBranchLength(childEvent) / 2;

	return event_construct2(NULL, branchLength, parentEvent, childEvent,
			eventTree);
}

/*
 * Creates the interpolation of a top node
 */
static Cap *buildFaces_interpolateTopNode(Face * face, int32_t topIndex) {
	Cap *topNode = face_getTopNode(face, topIndex);
	Cap *derivedEdgeDestination = face_getDerivedDestination(face, topIndex);
	int32_t bottomNodeIndex;
	int32_t bottomNodeNumber = face_getBottomNodeNumber(face, topIndex);
	Cap *derivedEdgeBottomNode = NULL;
	Event *topEvent, *bottomEvent, *interpolatedEvent;

	// If no derived edge
	if (derivedEdgeDestination == NULL)
		return NULL;

	// Search for bottom node which generated the derived edge
	for (bottomNodeIndex = 0; bottomNodeIndex < bottomNodeNumber; bottomNodeIndex++) {
		if (face_getDerivedDestinationAtIndex(face, topIndex, bottomNodeIndex)) {
			derivedEdgeBottomNode = face_getBottomNode(face, topIndex,
					bottomNodeIndex);
			break;
		}
	}

#ifdef BEN_DEBUG
	assert(derivedEdgeBottomNode);
#endif

	// Go to the appropriate descent edge
	while (cap_getParent(derivedEdgeBottomNode) != topNode)
		derivedEdgeBottomNode = cap_getParent(derivedEdgeBottomNode);

	// Event interpolation
	topEvent = cap_getEvent(topNode);
	bottomEvent = cap_getEvent(derivedEdgeBottomNode);
	while (event_getParent(bottomEvent) != topEvent)
		bottomEvent = event_getParent(bottomEvent);
	interpolatedEvent = buildFaces_interpolateEvents(topEvent, bottomEvent);

	// Cap interpolation
	return buildFaces_interpolateCaps(topNode, derivedEdgeBottomNode,
			interpolatedEvent);
}

/*
 * Produces array of interpolations for the top nodes of a face
 */
static Cap **buildFaces_interpolateTopNodes(Face * face) {
	Cap **interpolations = callocLocal(face_getCardinal(face), sizeof(Cap *));
	uint32_t topIndex;

	for (topIndex = 0; topIndex < face_getCardinal(face); topIndex++)
		interpolations[topIndex]
				= buildFaces_interpolateTopNode(face, topIndex);

	return interpolations;
}

/*
 * Connects an interpolated node within a face
 */
static void buildFaces_connectInterpolatedNode(Cap ** interpolations,
		int32_t nodeIndex, Face * face) {
	Cap *node = interpolations[nodeIndex];
	Cap *topNode = face_getTopNode(face, nodeIndex);
	Cap *adjacentNode = NULL;
	int32_t adjacentIndex;

	if (cap_getAdjacency(topNode))
		adjacentNode = cap_getPositiveOrientation(cap_getAdjacency(topNode));

	// If not interpolated or previously connected or top node disconnected
	if (node == NULL || cap_getAdjacency(node) || adjacentNode == NULL)
		return;

	// Look for index of adjacent node
	for (adjacentIndex = 0; adjacentIndex < face_getCardinal(face); adjacentIndex++)
		if (face_getTopNode(face, adjacentIndex) == adjacentNode)
			break;

#ifdef BEN_DEBUG
	// What if adjacent node is not a top node??
	if (adjacentIndex == face_getCardinal(face))
		abort();
#endif

	// Tie the knot
	cap_makeAdjacent(node, interpolations[adjacentIndex]);
}

/*
 * Connects the interpolated nodes of a face
 */
static void buildFaces_connectInterpolatedNodes(Cap ** interpolations,
		Face * face) {
	int32_t topIndex;

	for (topIndex = 0; topIndex < face_getCardinal(face); topIndex++)
		buildFaces_connectInterpolatedNode(interpolations, topIndex, face);
}

/*
 * Fill in interpolated face for a given top node
 */
static void buildFaces_fillInterpolatedFace(Face * face, Cap ** interpolations,
		Face * interpolatedFace, int32_t nodeCount, int32_t nodeIndex) {
	int32_t index, bottomNodeIndex;
	Cap * derivedDestination = face_getDerivedDestination(face, nodeIndex);

	// Find interpolation of derived destination
	for (index = 0; index < face_getCardinal(face); index++)
		if (face_getTopNode(face, index) == derivedDestination)
			break;

	// Fill face
	face_setTopNode(interpolatedFace, nodeCount, interpolations[nodeIndex]);
	face_setBottomNodeNumber(interpolatedFace, nodeCount, 1);
	face_setDerivedDestination(interpolatedFace, nodeCount, 0,
			interpolations[index]);

	// Search for bottom node which generated the derived edge
	for (bottomNodeIndex = 0; bottomNodeIndex < face_getBottomNodeNumber(face,
			nodeIndex); bottomNodeIndex++) {
		if (face_getDerivedDestinationAtIndex(face, nodeCount, bottomNodeIndex)) {
			face_addBottomNode(interpolatedFace, nodeCount,
					face_getBottomNode(face, nodeCount, bottomNodeIndex));
			break;
		}
	}
}

/*
 * Create face from interpolated nodes
 */
static void buildFaces_createInterpolatedFace(Face * face,
		Cap ** interpolations, Net * net) {
	Face *interpolatedFace = face_construct(net);
	int32_t nodeIndex;
	int32_t nodeCount = 0;

	// Count interpolated top nodes in face
	for (nodeIndex = 0; nodeIndex < face_getCardinal(face); nodeIndex++)
		if (interpolations[nodeIndex] != NULL)
			nodeCount++;

	// Initialize face
	face_allocateSpace(face, nodeCount);

	// Project face info onto interpolated face
	nodeCount = 0;
	for (nodeIndex = 0; nodeIndex < face_getCardinal(face); nodeIndex++)
		if (interpolations[nodeIndex] != NULL)
			buildFaces_fillInterpolatedFace(face, interpolations,
					interpolatedFace, nodeCount++, nodeIndex);
}

/*
 * Isolates into regular and trivial faces
 */
void buildFaces_isolate(Face * face, Net * net) {
	Cap **interpolations;

	// If uncessary
	if (face_isRegular(face))
		return;

	// Interpolate top nodes
	interpolations = buildFaces_interpolateTopNodes(face);

	// Connect intepolated nodes
	buildFaces_connectInterpolatedNodes(interpolations, face);

	// Create rearrangment face
	buildFaces_createInterpolatedFace(face, interpolations, net);

	// Cleaning up
	face_destruct(face);
	free(interpolations);

#ifdef BEN_DEBUG
	assert(face_isRegular(face));
#endif
}

/*
 * Isolates all the faces in the net
 */
void buildFaces_isolateFaces(Net * net) {
	Net_FaceIterator *iter = net_getFaceIterator(net);
	Face *face;

	logInfo("Isolating faces\n");

	while ((face = net_getNextFace(iter)))
		buildFaces_isolate(face, net);

	net_destructFaceIterator(iter);
}

/*
 * Builds a stub node and creates at the same time an ad hoc end
 */
static Cap * buildFaces_constructStub(Cap * adjacentCap) {
	//Construct the new stub and the new cap..
	Net * net = end_getNet(cap_getEnd(adjacentCap));
	End *newFreeStubEnd = createNewFreeStubEnd(net);
	Cap *cap = cap_construct(newFreeStubEnd, cap_getEvent(adjacentCap));

	//Now set the group of the new stub end (they should be in the same group)
	End *adjacentCapEnd = cap_getEnd(adjacentCap);
	Group *group = end_getGroup(adjacentCapEnd);
	end_setGroup(newFreeStubEnd, group);

	//Now make adjacent
	cap_makeAdjacent(cap, adjacentCap);

	return cap;

	/*
	 The 0 argument to the end constructor is a bool saying the stub end is 'free',
	 ie. not necessarily inherited from the parent (though it can be), and not part
	 of the reference structure. I.e. some stub ends are 'attached' - opposite of free,
	 representing the case where we know what happened to the other end of the stub
	 (i.e. if it is a block end at a higher level or the defined end of a sequence defined
	 at the top level).
	 */
}

/*
 * Engineers a node so that a regular face has an even number of
 * top/bottom node pairs
 */
static void buildFaces_engineerCaps(Face * face, Net * net) {
	int32_t index;
	Cap *nonAdjacent = NULL;
	int32_t nonDerived = -1;
	Cap *X, *Xprime;
	Cap *parent;

	// Look for terminal nodes in face according to their nature
	for (index = 0; index < face_getCardinal(face); index++) {
		if (!cap_getAdjacency(face_getTopNode(face, index))) {
			nonAdjacent = face_getTopNode(face, index);
			if (nonDerived > -1)
				break;
		} else if (face_getDerivedDestination(face, index) == NULL) {
			nonDerived = index;
			if (nonAdjacent)
				break;
		}
	}

	// Engineer appropriate nodes
	X = buildFaces_constructStub(nonAdjacent);

	Xprime = cap_construct(cap_getEnd(X), cap_getEvent(face_getBottomNode(face,
			nonDerived, 0)));
	cap_makeParentAndChild(X, Xprime);
	cap_makeAdjacent(Xprime, face_getBottomNode(face, nonDerived, 0));

	if (cap_getParent(nonAdjacent)) {
		parent = cap_construct(cap_getEnd(X), cap_getEvent(cap_getParent(
				nonAdjacent)));
		cap_makeParentAndChild(parent, X);
	}

	// Correct face:
	face_engineerArtificialNodes(face, X, Xprime, nonDerived);
}

/*
 * Connects the ends of an even length regular path to from
 * a canonical cycle
 */
static void buildFaces_close(Face * face) {
	int32_t index;
	Cap *nonAdjacent1 = NULL;
	Cap *nonAdjacent2 = NULL;
	int32_t nonDerived1 = -1;
	int32_t nonDerived2 = -1;

	// Look for terminal nodes in face and determine their nature
	for (index = 0; index < face_getCardinal(face); index++) {
		if (!cap_getAdjacency(face_getTopNode(face, index))) {
			if (nonAdjacent1 == NULL)
				nonAdjacent1 = face_getTopNode(face, index);
			else {
#ifndef BEN_DEBUG
				nonAdjacent2 = face_getTopNode(face, index);
				break;
#else
				// What if more than two nodes with no adjacency?
				if (nonAdjacent2)
					abort();
				else
					nonAdjacent2 = face_getTopNode(face, index);
#endif
			}
		} else if (face_getDerivedDestination(face, index) == NULL) {
			if (nonDerived1 == -1)
				nonDerived1 = index;
			else {
#ifndef BEN_DEBUG
				nonDerived2 = index;
				break;
#else
				// What if more than two nodes with no derived edge?
				if (nonDerived2 > -1)
					abort();
				else
					nonDerived2 = index;
#endif
			}
		}
	}

#ifdef BEN_DEBUG
	// What if some nodes with no adjacency and others with no derived edge?
	if (nonAdjacent1 && nonDerived1 > 1)
		abort();
#endif

	// If the terminal nodes lack an adjacency edge
	if (nonAdjacent2)
		cap_makeAdjacent(nonAdjacent1, nonAdjacent2);
	else if (nonDerived2 > -1)
		cap_makeAdjacent(face_getBottomNode(face, nonDerived1, 0),
				face_getBottomNode(face, nonDerived2, 0));
}

/*
 * Canonizes face into a regular cycle
 */
void buildFaces_canonize(Face * face, Net * net) {
	if (face_getCardinal(face) % 2 == 0)
		buildFaces_engineerCaps(face, net);

#ifdef BEN_DEBUG
	assert(face_isCanonical(face));
#endif 

	buildFaces_close(face);
}

/*
 * Canonizes all the faces in the net
 */
void buildFaces_canonizeFaces(Net * net) {
	Net_FaceIterator *iter = net_getFaceIterator(net);
	Face *face;

	logInfo("Canonizing faces\n");

	while ((face = net_getNextFace(iter)))
		if (!face_isCanonical(face))
			buildFaces_canonize(face, net);

	net_destructFaceIterator(iter);
}

/*
 * The works: create, regularize and canonize faces in net
 */
void buildFaces_buildAndProcessFaces(Net * net) {
	buildFaces_constructFaces(net);
	buildFaces_simplifyFaces(net);
	buildFaces_isolateFaces(net);
	buildFaces_canonizeFaces(net);
}

//Todo - fine me a home.
End *createNewFreeStubEnd(Net *net) {
	End *newFreeStubEnd = end_construct(0, net);
	assert(net_getGroupNumber(net) == 1);
	end_setGroup(newFreeStubEnd, net_getFirstGroup(net));
	return newFreeStubEnd;
}

