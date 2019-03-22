set.seed(1);
nodeId = tree_sample$node;
parent = tree_sample$parent;
d2p = tree_sample$distance;
nLeaves = 4;
nNodes = 8;
leaves = 1:nLeaves;
father = 0;
distanceToFather = 0;
lambda0 = 0.2;
lambda_dup = 0.1;
lambfa_trans = 0.01;
lambda_loss = 0.05;
set = list();
for (i in 1:nLeaves){set[[i]] = c(i)};
lable = rep(0, nNodes);
root = 7;

findFather <- function(nodeId){
  return(parent[nodeId])
}
findChild <- function(fatherId){
  child = which (parent == fatherId)
  return(child)
}
findDistance <- function(nodeId){
  return(d2p[nodeId])
}
mark = rep(0,nNodes)
mark[1:nLeaves] = 1

coal <- function(node, mark, distance, set1){
  #print(distance)
  if (mark <= 1) return(list("mark" = mark,"set" = set1, "distance" = distance))
  else {
    lambda = mark * lambda0
    distance_fake = rexp(1, rate = lambda)
    if (distance < distance_fake) return(list("mark" = mark,"set" = set1, "distance" = distance))
    else {
      # change set
      if (length(set1[[node]]) >= 2) {
        print(paste0("initial ", "node ", node, " : ", paste0(set1[[node]], collapse = ", ")));
        couple = sample(set1[[node]], size = 2, replace = FALSE)
        for (i in 1:length(couple)){
          set1[[node]] = set1[[node]][!(set1[[node]] == couple[i])]
        }
        new_entry = paste0(couple[1],couple[2])
        set1[[node]] = append(new_entry, set1[[node]], 0)
        # print out node info
        print(paste0("coalescent at ", "node ", node, " : ",
                     paste0(set1[[node]], collapse = ", "), ", distance = ", distance_fake));
      }
      else{
        return(list("mark" = mark,"set" = set1, "distance" = distance))
      }
      distance = distance - distance_fake
      coal(node, mark-1, distance, set1)
    }
  }
}

new_leaves = NULL

while(1){
  for (leaf in leaves) {
    if (leaf == root) {
      tempList = coal(root, mark[root], 10000, set)
      # temp = set[[leaf]][1]
      # if (length(set[[leaf]]) >= 2) {
      #   for (i in 2:length(set[[leaf]])) {
      #     temp = paste0(temp, set[[leaf]][i])
      #   }
      #   set[[leaf]] = temp
      # }
      break
    }
    else{
      father = findFather(leaf)
      child = findChild(father)
      if (lable[child[1]] == 1) next
      d1 = findDistance(child[1])
      d2 = findDistance(child[2])
      lable[leaf] = 1;
      if (mark[child[1]] != 0 && mark[child[2]] != 0){
        tempList1 = coal(child[1], mark[child[1]], d1, set)
        tempList2 = coal(child[2], mark[child[2]], d2, set)
        mark[father] = tempList1$mark
        set[[child[1]]] = tempList1$set[[child[1]]]
        mark[father] = mark[father] + tempList2$mark
        set[[child[2]]] = tempList2$set[[child[2]]]
        set[[father]] = union(set[[child[1]]],set[[child[2]]])
        if (length(new_leaves) > 0) {
          for (i in 1:length(new_leaves)) {
            if (new_leaves[i] == child[1]) {
              new_leaves = new_leaves[!which(new_leaves == child[1])]
            }
            if (new_leaves[i] == child[2]) {
              new_leaves = new_leaves[!which(new_leaves == child[2])]
            }
          }
        }
        new_leaves = append(father, new_leaves, 0)
      }
      else {
        new_leaves = append(leaf, new_leaves, 0)
      }
    }
  }
  if (leaf == root) break
  new_leaves = new_leaves[!duplicated(new_leaves)]
  leaves = new_leaves
  lable = rep(0, nNodes)
  new_leaves = NULL
}