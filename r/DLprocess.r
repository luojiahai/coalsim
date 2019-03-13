geneTree <- read.csv("./geneTree.csv")

node = geneTree$node;
child1 = geneTree$child1;
child2 = geneTree$child2;
d_child1 = geneTree$d_child1;
d_child2 = geneTree$d_child2;
lambda_d = 1
lambda_l = 0.3
root = 7

scan <- function(node, distance){ #print distance between the event node and its child
  d_dup = rexp(1, lambda_d)
  d_loss = rexp(1, lambda_l)
  if (d_dup < d_loss && d_dup < distance) {
    # duplication first, continue
    # adding the duplicated subtree using coalesence
    print(paste0("duplication at node ", node, " with distance ", distance - d_dup))
    scan(node, distance - d_dup)
  }
  else if (d_loss <= d_dup && d_loss < distance ) {
    # loss first, break
    print(paste0("loss at node ", node, " with distance ", distance - d_loss))
  }
  else { # when both greater than distance
    # nothing happened, continue
    print(paste0("nothing happened at node ", node))
    kid1 = child1[node]
    kid2 = child2[node]
    distance_to_kid1 = d_child1[node]
    distance_to_kid2 = d_child2[node]
    if (kid1 != -1) {
      scan(kid1, distance_to_kid1)
      scan(kid2, distance_to_kid2)
    }
    else{
      print(paste0("reach the end of node ", node))
    }
  }
}

scan(root, 0)