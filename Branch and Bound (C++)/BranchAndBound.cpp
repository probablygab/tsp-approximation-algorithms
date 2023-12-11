#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <cmath>
#include <queue>
#include <bitset>
#include <stack>
#include <chrono>
#include <climits>

// Needed for python bindings
#include <pybind11/pybind11.h>

// Max number of nodes in a graph (needed for bitmask)
#define MAX_NODES 64

// Memory optimized struct (12 + MAX_NODES/8) total bytes
// O(1) everything, with 64 nodes we have 20 bytes per node
struct Node {
    int bound;
    int level;
    int lastNode;
    std::bitset<MAX_NODES> bitmask;
    
    Node(int bnd, int lvl, std::bitset<MAX_NODES> msk, int last) {
        bound = bnd;
        level = lvl;
        lastNode = last;

        bitmask = msk;
        bitmask.set(last);
    };
};

struct SmallestEdges {
    int smallestWeight       = INT_MAX;
    int secondSmallestWeight = INT_MAX;
};

// Literally just used so we can build a graph
struct Point {
    double x;
    double y;
};

typedef std::vector<std::vector<int>> Graph;
typedef std::vector<Point> vecp;
typedef std::vector<SmallestEdges> vecedge;

void buildGraph(Graph &graph, vecedge &smallEdges, const char* filename) {
    FILE* file = fopen(filename, "rb");

    static char buffer[200];
    int size;
    int lineIdx = 0;
    int read = 0;
    int trash;

    // Vector for points
    vecp nodes;

    while (true) {
        fgets(buffer, 196, file);
        lineIdx++;

        // Get dimension
        if (lineIdx == 4) {
            sscanf(&buffer[10], "%d", &size);

            graph.resize(size);

            for (std::vector<int> &v : graph)
                v.resize(size);

            nodes.resize(size);
        }

        // Start to get nodes
        if (lineIdx > 6) {
            sscanf(buffer, "%d %lf %lf\n", &trash, &nodes[read].x, &nodes[read].y);
            read++;
        }
        
        if (read == size)
            break;
    }

    // Create graph edges
    // Already calculate two smallest edges for each node
    smallEdges.resize(size);

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            // Skip edge to self
            if (i == j)
                continue;

            double xd = nodes[i].x - nodes[j].x;
            double yd = nodes[i].y - nodes[j].y;

            int dist = static_cast<int>(std::round(std::sqrt(xd*xd + yd*yd)));

            // Check for small edges
            if (dist < smallEdges[i].smallestWeight) {
                smallEdges[i].secondSmallestWeight = smallEdges[i].smallestWeight;
                smallEdges[i].smallestWeight = dist;
            } else if (dist < smallEdges[i].secondSmallestWeight) {
                smallEdges[i].secondSmallestWeight = dist;
            }

            graph[i][j] = dist;
        }

    fclose(file);
}

int getBaseBound(const vecedge &smallEdges) {
    int bound = 0;

    for (const SmallestEdges &edges : smallEdges)
        bound += edges.smallestWeight + edges.secondSmallestWeight;

    return bound;
}

inline int getBound(const Node &base, const Graph &graph, const vecedge &smallEdges, const int k) {
    // Compute bound in O(1) thanks to precomputing stuff yay :)
    int bound = base.bound;
    int last = base.lastNode;

    // Add new edge weights
    bound += 2 * graph[last][k];

    // Remove smallest edge from last (second is already taken)
    if (last != 0)
        bound -= smallEdges[last].smallestWeight;
    else
        bound -= smallEdges[last].secondSmallestWeight;

    // Remove 2nd smallest edge from current (just one edge is taken)
    if (k != 0)
        bound -= smallEdges[k].secondSmallestWeight;
    else
        bound -= smallEdges[k].smallestWeight;
    
    return bound;
}

inline int convBound(int bound) {
    return static_cast<int>(std::ceil((double) bound / 2.0));
}

int branchAndBound(int best, const char* filename, int timeoutSeconds) {
    // Get graph edges and precompute smallest edges and bound
    Graph graph;
    vecedge smallEdges;
    buildGraph(graph, smallEdges, filename);
    int base = getBaseBound(smallEdges);

    // Init timeout (after building graph)
    auto endTime = std::chrono::system_clock::now() + std::chrono::seconds(timeoutSeconds);   

    // Init pq using a lambda for comparations
    auto cmp = [](const Node &left, const Node &right) { return left.bound < right.bound; };
    std::priority_queue<Node, std::vector<Node>, decltype(cmp)> pq(cmp);

    // Init stack to simulate iterative DFS
    Node root = Node(base, 1, 0, 0);

    std::stack<Node> st;
    st.push(root);

    // Main loop
    while (!st.empty()) {
        // Check timeout
        if (std::chrono::system_clock::now() > endTime)
            return best;

        Node node = st.top();
        st.pop();

        // Get cost from bound
        // if all nodes are covered the bound is exactly the cost of the solution
        int cost = convBound(node.bound);
        
        if (cost > best)
            continue;

        // Leaf node, possible solution
        if (node.level > graph.size()) {
            if (cost < best)
                best = cost;
        // Not a leaf node, explore solutions if possible
        } else if (cost < best) {
            if (node.level < graph.size()) {
               
                for (int k = 1; k < graph.size(); k++) {
                    // If node k not in solution
                    if (!node.bitmask.test(k)) {
                        int lowerBound = getBound(node, graph, smallEdges, k);

                        // Bound is promising, add node to queue
                        if (convBound(lowerBound) < best) {
                            Node newNode = Node(lowerBound, 
                                           node.level + 1, 
                                           node.bitmask,
                                           k);
                            pq.push(newNode);
                        }
                    }
                }

                // Add nodes from pq to stack
                // thus nodes with lower bounds will be unstacked first
                while (!pq.empty()) {
                    st.push(pq.top());
                    pq.pop();
                }
            // All levels covered, add return edge
            } else {
                // Return edge for first node (0)
                int lowerBound = getBound(node, graph, smallEdges, 0);

                if (convBound(lowerBound) < best) {
                    Node newNode = Node(lowerBound, 
                                    node.level + 1, 
                                    node.bitmask,
                                    0);
                    st.push(newNode);
                }
            }
        }
    }

    return best;
}

PYBIND11_MODULE(bnb, m) {
    m.def("branchAndBound",
          &branchAndBound, 
          pybind11::arg("best"),
          pybind11::arg("filename"),
          pybind11::arg("timeoutSeconds"),
          pybind11::return_value_policy::automatic,
          "Branch and Bound for TSP");
}

// int main(int argc, char* argv[]) {
//     if (argc < 4) {
//         printf("Usage: %s [problem file] [initial best] [timeout in seconds]\n", argv[0]);
//         exit(1);
//     }

//     int sol = branchAndBound(atoi(argv[2]), argv[1], atoi(argv[3]));
//     printf("Best: %d\n", sol);

//     return 0;
// }