library(randomForest)

library(sdmpredictors)
library(dismo)
library(biomod2)
library(ecospat)
library(caret)
library(pROC)


#devtools::install_github("sigopt/SigOptR")

layers <- c("BO_sstrange", "BO_sstmax", "BO_parmean","BO_salinity")
env <- load_layers(layers, equalarea = FALSE)

gf <- function(f, native = TRUE) {
  if(native) {
    d <- read.csv(file.path("data", "undaria_native_Europe_spthin100_coastal", f))
  } else {
    d <- read.csv(file.path("data", "undaria_all_Europe_spthin100_coastal", f))
  }
  colnames(d) <- c("x", "y")
  d[,1:2]
}

create_model <- function(occ_train, occ_test, bg_train, bg_test, ...) {
  train <- rbind(occ_train, bg_train)
  train_y <- c(rep(1, nrow(occ_train)), rep(0, nrow(bg_train)))

  test <- rbind(occ_test, bg_test)
  test_y <- c(rep(1, nrow(occ_test)), rep(0, nrow(bg_test)))

  #model <- randomForest::randomForest(train, as.factor(train_y), mtry=2, nodesize = 5)
  #model <- randomForest::randomForest(train, as.factor(train_y))
  model <- randomForest::randomForest(train, as.factor(train_y), ...)

  #training_fit <- predict(model, type="prob") ## Out-of-bag prediction probablities
  #pROC::auc(train_y, training_fit[,2])
  #dismo::evaluate(occ_test, bg_test, model, type="prob")
  test_fit <- predict(model, test, as.factor(test_y), type="prob")
  test_auc <- pROC::auc(test_y, test_fit[,2])


  stats <- lapply(1:(ceiling(1000*max(test_fit))), function(tr) {
    misc <- table(as.integer(test_fit[,2] >= tr/1000), test_y)
    data.frame(treshold=tr/1000,
               kappa=biomod2::calculate.stat(misc, stat='KAPPA'),
               tss=biomod2::calculate.stat(misc, stat='TSS'))
  })
  stats <- do.call(rbind, stats)
  kappa <- stats[which.max(stats$kappa), c("kappa", "treshold")]
  colnames(kappa) <- c("stat", "treshold")
  tss <- stats[which.max(stats$tss), c("tss", "treshold")]
  colnames(tss) <- c("stat", "treshold")
  evaluation <- rbind(cbind(metric="kappa", kappa),
                      cbind(metric="tss", tss),
                      data.frame(metric="auc", stat=test_auc, treshold=NA))
  ## TODO save models
  evaluation
}
rep_many_small <- function(occ_train, occ_test, bg_train, bg_test, ..., nrep = 5, smallsize = 0) {

}

rep_models <- function(..., nrep = 10) {
  evals <- data.frame(row.names = c("kappa", "tss", "auc"))

  for(i in 1:nrep) {
    r <- create_model(...)
    evals <- cbind(evals, rep=r$stat*100)
  }
  ee <- as.data.frame(t(evals))
  p <- ggplot(ee, aes(x=auc)) + geom_density() + geom_point(y=0)
  print(p)
  e <- cbind(mean = rowMeans(evals),
             sd = apply(evals, 1, sd),
             min = apply(evals, 1, min),
             max = apply(evals, 1, max))
  print(e)
  e
}

occ_train_xy <- gf("occurrences_train.csv", native = TRUE)
occ_test_xy <- gf("occurrences_test.csv", native = TRUE)
bg_test_xy <- gf("background_test.csv")
bg_train_xy <- read.csv("data/coastal_background_train.csv")

# XY null model
xy_eval <- create_model(occ_train_xy, occ_test_xy, bg_train_xy, bg_test_xy)

occ_train <- extract(env, occ_train_xy)
occ_test <- extract(env, occ_test_xy)
bg_train <- extract(env, bg_train_xy)
bg_test <- extract(env, bg_test_xy)

native_eval <- create_model(occ_train, occ_test, bg_train, bg_test)
native_eval_replace_false <- create_model(occ_train, occ_test, bg_train, bg_test, replace = FALSE)
native_eval_m2_n5 <- create_model(occ_train, occ_test, bg_train, bg_test, mtry=2, nodesize = 5)

create_model(occ_train, occ_test, bg_train, bg_test, classwt=c(0.5,0.5))
create_model(occ_train, occ_test, bg_train, bg_test, maxnodes = 100)
create_model(occ_train, occ_test, bg_train, bg_test, maxnodes = 100, ntree = 2000)
create_model(occ_train, occ_test, bg_train, bg_test, maxnodes = 500)
rep_models(occ_train, occ_test, bg_train, bg_test, maxnodes = 500)
create_model(occ_train, occ_test, bg_train, bg_test, maxnodes = 500, ntree = 1000)
create_model(occ_train, occ_test, bg_train, bg_test, maxnodes = 1000)
create_model(occ_train, occ_test, bg_train, bg_test, ntree = 1000)


create_model(occ_train, occ_test, bg_train, bg_test, sampsize = c(30, 30))
create_model(occ_train, occ_test, bg_train, bg_test, sampsize = c(70, NROW(occ_train)))

create_model(occ_train, occ_test, bg_train, bg_test, sampsize = c(NROW(occ_train), NROW(occ_train)))
sampsize_nrow <- rep_models(occ_train, occ_test, bg_train, bg_test, sampsize = c(NROW(occ_train), NROW(occ_train)))
sampsize_nrow_nodesize5 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                        nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_ntree2000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                ntree=2000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_ntree3000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                      ntree=3000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_ntree4000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                ntree=4000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_ntree5000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                ntree=5000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

## best ???
sampsize_nrow_nodesize5_ntree10000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                ntree=10000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize3_ntree10000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                ntree=10000, nodesize = 3, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize10 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                       nodesize = 10, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize10_ntree1000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                       ntree = 1000, nodesize = 10, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize50 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                       nodesize = 50, sampsize = c(NROW(occ_train), NROW(occ_train)))


rep_models(occ_train, occ_test, bg_train, bg_test,
           nodesize = 50, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize50_ntree2000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                       ntree = 2000, nodesize = 50, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_mtry1 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                            mtry=1, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_mtry2 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                            mtry=2, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_mtry3 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                            mtry=3, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

sampsize_nrow_nodesize5_mtry4 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                            mtry=4, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))

native_eval_sampsize35_ntree2000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                               sampsize = c(35, 35), ntree = 2000)

create_model(occ_train, occ_test, bg_train, bg_test, sampsize = c(300, NROW(occ_train)))


rep_models(occ_train, occ_test, bg_train, bg_test,
           nodesize = 2,
           maxnodes = 5200, ntree = 10, nrep=1)


## best ???
sampsize_nrow_nodesize5_ntree10000 <- rep_models(occ_train, occ_test, bg_train, bg_test,
                                                 ntree=10000, nodesize = 5, sampsize = c(NROW(occ_train), NROW(occ_train)))



experiment <- SigOptR::create_experiment(list(
  name="rf4sdm undaria native VS europe",
  parameters=list(
    list(name="sampsize", type="int", bounds=list(min=30, max=NROW(occ_train))),
    list(name="ntree", type="int", bounds=list(min=1, max=10)),
    list(name="nodesize", type="int", bounds=list(min=1, max=20)),
    list(name="replace", type="int", bounds=list(min=0, max=1)),
    list(name="maxnodes_perc", type="int", bounds=list(min=70, max=100)) ## maxnodes percentage
  )
))
# experiment_id = 3977
experiment_id = experiment$id # 3977
results <- list()

for(i in 1:50) {
  suggestion <- SigOptR::create_suggestion(experiment_id)
  p <- suggestion$assignments
  tryCatch({
    res <- rep_models(occ_train, occ_test, bg_train, bg_test,
                      nodesize = p$nodesize, ntree = p$ntree*1000,
                      replace = as.logical(p$replace), sampsize = c(p$sampsize,p$sampsize),
                      maxnodes = trunc((trunc(p$sampsize/p$nodesize)+1) * (p$maxnodes_prec/100)))

    results[[suggestion$id]] <- list(p, res)

    #
    SigOptR::create_observation(experiment_id,
                                list(suggestion=suggestion$id,
                                     value=res["auc","mean"],
                                     value_stddev=res["auc","sd"]))
  }, error = function(e) {
    print("FAILED")
    SigOptR::create_observation(experiment_id, list(suggestion=suggestion$id, failed=TRUE))
  })
}
