#include "clique.h"
int main() {
    std::vector<double> w = {2, 2, 4};
    std::vector<double> e = {5, 5, 5};
    CliqueCenter graph(w, e, 2);
    graph.printCenter();
}