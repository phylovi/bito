// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
//

#include "subsplit_dag_storage.hpp"
#include "nni_operation.hpp"

class TPComputationEngine {
 public:
  //
  virtual void Initialize();
  //
  virtual double GetTopTreeScoreWithProposedNNI(const NNIOperation &post_nni,
                                                const NNIOperation &pre_nni);
  //
  virtual double GetTopTreeScoreWithEdge(const EdgeId edge_id);
  //
  virtual void Compute();
  //
  virtual void UpdateAfterAddNodePair();
}
