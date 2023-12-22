Code needed to run ARTEMIS and DELFI in the Lung detection

Each model subfolder contains subfolders Cohorts and Control. The Cohorts folder contains the machine learning features, the Control folder contains code to train (and if applicable, test) models.

In the Control folder are trained by calling LOO.sh which calls Ensemble_ARTEMIS1_delfi_LOO.r, and generates cross-validated results. These may be summarized with performance.r. In the Lung_Detection model, Lock.sh which calls Train_Lock.r is also provided to externally validate the model on a test-set.  

The code provided trains several different model architectures using repeat landscape (ARTEMIS) and fragmentation (DELFI) features -- the final models chosen for publication are described in the methods section of the published paper and documented in the code used to generate relevant figures and tables. The final model scores for each sample used are available in the supplementary methods of the paper.

 

