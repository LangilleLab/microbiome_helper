#####################################################################
# Author: Jacob Nearing;
# Date April 10th 2019
# Title: Script to run the main RF pipeline
#####################################################################
#
###### IMPORT DEPENDECIES ###########################################

deps = c("randomForest", "pROC", "caret", "DMwR", "doMC", "tidyverse", "PRROC", "MLmetrics", "reshape2")
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = "http://cran.us.r-project.org")
  }
  library(dep, character.only = TRUE)
}

######################## GLOBAL VARIABLES ###########################
CORES_TO_USE <- 1

COLORS <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", 
            "#008941", "#006FA6", "#A30059",
            "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", 
            "#004D43", "#8FB0FF", "#997D87",
            "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", 
            "#3B5DFF", "#4A3B53", "#FF2F80",
            "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", 
            "#FF90C9", "#B903AA", "#D16100",
            "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", 
            "#0AA6D8", "#013349", "#00846F",
            "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", 
            "#C0B9B2", "#C2FF99", "#001E09",
            "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", 
            "#B77B68", "#7A87A1", "#788D66",
            "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", 
            "#456648", "#0086ED", "#886F4C",
            "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", 
            "#A3C8C9", "#FF913F", "#938A81",
            "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", 
            "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", 
            "#772600", "#D790FF", "#9B9700",
            "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", 
            "#99ADC0", "#3A2465", "#922329",
            "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
            "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
            "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
            "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
            "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
            "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
            "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
            "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
            "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
            "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
            "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
            "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
            "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
            "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
            "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
            "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
            "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
            "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B", "#1E2324", "#DEC9B2", "#9D4948",
            "#85ABB4", "#342142", "#D09685", "#A4ACAC", "#00FFFF", "#AE9C86", "#742A33", "#0E72C5",
            "#AFD8EC", "#C064B9", "#91028C", "#FEEDBF", "#FFB789", "#9CB8E4", "#AFFFD1", "#2A364C",
            "#4F4A43", "#647095", "#34BBFF", "#807781", "#920003", "#B3A5A7", "#018615", "#F1FFC8",
            "#976F5C", "#FF3BC1", "#FF5F6B", "#077D84", "#F56D93", "#5771DA", "#4E1E2A", "#830055",
            "#02D346", "#BE452D", "#00905E", "#BE0028", "#6E96E3", "#007699", "#FEC96D", "#9C6A7D",
            "#3FA1B8", "#893DE3", "#79B4D6", "#7FD4D9", "#6751BB", "#B28D2D", "#E27A05", "#DD9CB8",
            "#AABC7A", "#980034", "#561A02", "#8F7F00", "#635000", "#CD7DAE", "#8A5E2D", "#FFB3E1",
            "#6B6466", "#C6D300", "#0100E2", "#88EC69", "#8FCCBE", "#21001C", "#511F4D", "#E3F6E3",
            "#FF8EB1", "#6B4F29", "#A37F46", "#6A5950", "#1F2A1A", "#04784D", "#101835", "#E6E0D0",
            "#FF74FE", "#00A45F", "#8F5DF8", "#4B0059", "#412F23", "#D8939E", "#DB9D72", "#604143",
            "#B5BACE", "#989EB7", "#D2C4DB", "#A587AF", "#77D796", "#7F8C94", "#FF9B03", "#555196",
            "#31DDAE", "#74B671", "#802647", "#2A373F", "#014A68", "#696628", "#4C7B6D", "#002C27",
            "#7A4522", "#3B5859", "#E5D381", "#FFF3FF", "#679FA0", "#261300", "#2C5742", "#9131AF",
            "#AF5D88", "#C7706A", "#61AB1F", "#8CF2D4", "#C5D9B8", "#9FFFFB", "#BF45CC", "#493941",
            "#863B60", "#B90076", "#003177", "#C582D2", "#C1B394", "#602B70", "#887868", "#BABFB0",
            "#030012", "#D1ACFE", "#7FDEFE", "#4B5C71", "#A3A097", "#E66D53", "#637B5D", "#92BEA5",
            "#00F8B3", "#BEDDFF", "#3DB5A7", "#DD3248", "#B6E4DE", "#427745", "#598C5A", "#B94C59",
            "#8181D5", "#94888B", "#FED6BD", "#536D31", "#6EFF92", "#E4E8FF", "#20E200", "#FFD0F2",
            "#4C83A1", "#BD7322", "#915C4E", "#8C4787", "#025117", "#A2AA45", "#2D1B21", "#A9DDB0",
            "#FF4F78", "#528500", "#009A2E", "#17FCE4", "#71555A", "#525D82", "#00195A", "#967874",
            "#555558", "#0B212C", "#1E202B", "#EFBFC4", "#6F9755", "#6F7586", "#501D1D", "#372D00",
            "#741D16", "#5EB393", "#B5B400", "#DD4A38", "#363DFF", "#AD6552", "#6635AF", "#836BBA",
            "#98AA7F", "#464836", "#322C3E", "#7CB9BA", "#5B6965", "#707D3D", "#7A001D", "#6E4636",
            "#443A38", "#AE81FF", "#489079", "#897334", "#009087", "#DA713C", "#361618", "#FF6F01",
            "#006679", "#370E77", "#4B3A83", "#C9E2E6", "#C44170", "#FF4526", "#73BE54", "#C4DF72",
            "#ADFF60", "#00447D", "#DCCEC9", "#BD9479", "#656E5B", "#EC5200", "#FF6EC2", "#7A617E",
            "#DDAEA2", "#77837F", "#A53327", "#608EFF", "#B599D7", "#A50149", "#4E0025", "#C9B1A9",
            "#03919A", "#1B2A25", "#E500F1", "#982E0B", "#B67180", "#E05859", "#006039", "#578F9B",
            "#305230", "#CE934C", "#B3C2BE", "#C0BAC0", "#B506D3", "#170C10", "#4C534F", "#224451",
            "#3E4141", "#78726D", "#B6602B", "#200441", "#DDB588", "#497200", "#C5AAB6", "#033C61",
            "#71B2F5", "#A9E088", "#4979B0", "#A2C3DF", "#784149", "#2D2B17", "#3E0E2F", "#57344C",
            "#0091BE", "#E451D1", "#4B4B6A", "#5C011A", "#7C8060", "#FF9491", "#4C325D", "#005C8B",
            "#E5FDA4", "#68D1B6", "#032641", "#140023", "#8683A9", "#CFFF00", "#A72C3E", "#34475A",
            "#B1BB9A", "#B4A04F", "#8D918E", "#A168A6", "#813D3A", "#425218", "#DA8386", "#776133",
            "#563930", "#8498AE", "#90C1D3", "#B5666B", "#9B585E", "#856465", "#AD7C90", "#E2BC00",
            "#E3AAE0", "#B2C2FE", "#FD0039", "#009B75", "#FFF46D", "#E87EAC", "#DFE3E6", "#848590",
            "#AA9297", "#83A193", "#577977", "#3E7158", "#C64289", "#EA0072", "#C4A8CB", "#55C899",
            "#E78FCF", "#004547", "#F6E2E3", "#966716", "#378FDB", "#435E6A", "#DA0004", "#1B000F",
            "#5B9C8F", "#6E2B52", "#011115", "#E3E8C4", "#AE3B85", "#EA1CA9", "#FF9E6B", "#457D8B",
            "#92678B", "#00CDBB", "#9CCC04", "#002E38", "#96C57F", "#CFF6B4", "#492818", "#766E52",
            "#20370E", "#E3D19F", "#2E3C30", "#B2EACE", "#F3BDA4", "#A24E3D", "#976FD9", "#8C9FA8",
            "#7C2B73", "#4E5F37", "#5D5462", "#90956F", "#6AA776", "#DBCBF6", "#DA71FF", "#987C95",
            "#52323C", "#BB3C42", "#584D39", "#4FC15F", "#A2B9C1", "#79DB21", "#1D5958", "#BD744E",
            "#160B00", "#20221A", "#6B8295", "#00E0E4", "#102401", "#1B782A", "#DAA9B5", "#B0415D",
            "#859253", "#97A094", "#06E3C4", "#47688C", "#7C6755", "#075C00", "#7560D5", "#7D9F00",
            "#C36D96", "#4D913E", "#5F4276", "#FCE4C8", "#303052", "#4F381B", "#E5A532", "#706690",
            "#AA9A92", "#237363", "#73013E", "#FF9079", "#A79A74", "#029BDB", "#FF0169", "#C7D2E7",
            "#CA8869", "#80FFCD", "#BB1F69", "#90B0AB", "#7D74A9", "#FCC7DB", "#99375B", "#00AB4D",
            "#ABAED1", "#BE9D91", "#E6E5A7", "#332C22", "#DD587B", "#F5FFF7", "#5D3033", "#6D3800",
            "#FF0020", "#B57BB3", "#D7FFE6", "#C535A9", "#260009", "#6A8781", "#A8ABB4", "#D45262",
            "#794B61", "#4621B2", "#8DA4DB", "#C7C890", "#6FE9AD", "#A243A7", "#B2B081", "#181B00",
            "#286154", "#4CA43B", "#6A9573", "#A8441D", "#5C727B", "#738671", "#D0CFCB", "#897B77",
            "#1F3F22", "#4145A7", "#DA9894", "#A1757A", "#63243C", "#ADAAFF", "#00CDE2", "#DDBC62",
            "#698EB1", "#208462", "#00B7E0", "#614A44", "#9BBB57", "#7A5C54", "#857A50", "#766B7E",
            "#014833", "#FF8347", "#7A8EBA", "#274740", "#946444", "#EBD8E6", "#646241", "#373917",
            "#6AD450", "#81817B", "#D499E3", "#979440", "#011A12", "#526554", "#B5885C", "#A499A5",
            "#03AD89", "#B3008B", "#E3C4B5", "#96531F", "#867175", "#74569E", "#617D9F", "#E70452",
            "#067EAF", "#A697B6", "#B787A8", "#9CFF93", "#311D19", "#3A9459", "#6E746E", "#B0C5AE",
            "#84EDF7", "#ED3488", "#754C78", "#384644", "#C7847B", "#00B6C5", "#7FA670", "#C1AF9E",
            "#2A7FFF", "#72A58C", "#FFC07F", "#9DEBDD", "#D97C8E", "#7E7C93", "#62E674", "#B5639E",
            "#FFA861", "#C2A580", "#8D9C83", "#B70546", "#372B2E", "#0098FF", "#985975", "#20204C",
            "#FF6C60", "#445083", "#8502AA", "#72361F", "#9676A3", "#484449", "#CED6C2", "#3B164A",
            "#CCA763", "#2C7F77", "#02227B", "#A37E6F", "#CDE6DC", "#CDFFFB", "#BE811A", "#F77183",
            "#EDE6E2", "#CDC6B4", "#FFE09E", "#3A7271", "#FF7B59", "#4E4E01", "#4AC684", "#8BC891",
            "#BC8A96", "#CF6353", "#DCDE5C", "#5EAADD", "#F6A0AD", "#E269AA", "#A3DAE4", "#436E83",
            "#002E17", "#ECFBFF", "#A1C2B6", "#50003F", "#71695B", "#67C4BB", "#536EFF", "#5D5A48",
            "#890039", "#969381", "#371521", "#5E4665", "#AA62C3", "#8D6F81", "#2C6135", "#410601",
            "#564620", "#E69034", "#6DA6BD", "#E58E56", "#E3A68B", "#48B176", "#D27D67", "#B5B268",
            "#7F8427", "#FF84E6", "#435740", "#EAE408", "#F4F5FF", "#325800", "#4B6BA5", "#ADCEFF",
            "#9B8ACC", "#885138", "#5875C1", "#7E7311", "#FEA5CA", "#9F8B5B", "#A55B54", "#89006A",
            "#AF756F", "#2A2000", "#7499A1", "#FFB550", "#00011E", "#D1511C", "#688151", "#BC908A",
            "#78C8EB", "#8502FF", "#483D30", "#C42221", "#5EA7FF", "#785715", "#0CEA91", "#FFFAED",
            "#B3AF9D", "#3E3D52", "#5A9BC2", "#9C2F90", "#8D5700", "#ADD79C", "#00768B", "#337D00",
            "#C59700", "#3156DC", "#944575", "#ECFFDC", "#D24CB2", "#97703C", "#4C257F", "#9E0366",
            "#88FFEC", "#B56481", "#396D2B", "#56735F", "#988376", "#9BB195", "#A9795C", "#E4C5D3",
            "#9F4F67", "#1E2B39", "#664327", "#AFCE78", "#322EDF", "#86B487", "#C23000", "#ABE86B",
            "#96656D", "#250E35", "#A60019", "#0080CF", "#CAEFFF", "#323F61", "#A449DC", "#6A9D3B",
            "#FF5AE4", "#636A01", "#D16CDA", "#736060", "#FFBAAD", "#D369B4", "#FFDED6", "#6C6D74",
            "#927D5E", "#845D70", "#5B62C1", "#2F4A36", "#E45F35", "#FF3B53", "#AC84DD", "#762988",
            "#70EC98", "#408543", "#2C3533", "#2E182D", "#323925", "#19181B", "#2F2E2C", "#023C32",
            "#9B9EE2", "#58AFAD", "#5C424D", "#7AC5A6", "#685D75", "#B9BCBD", "#834357", "#1A7B42",
            "#2E57AA", "#E55199", "#316E47", "#CD00C5", "#6A004D", "#7FBBEC", "#F35691", "#D7C54A",
            "#62ACB7", "#CBA1BC", "#A28A9A", "#6C3F3B", "#FFE47D", "#DCBAE3", "#5F816D", "#3A404A",
            "#7DBF32", "#E6ECDC", "#852C19", "#285366", "#B8CB9C", "#0E0D00", "#4B5D56", "#6B543F",
            "#E27172", "#0568EC", "#2EB500", "#D21656", "#EFAFFF", "#682021", "#2D2011", "#DA4CFF",
            "#70968E", "#FF7B7D", "#4A1930", "#E8C282", "#E7DBBC", "#A68486", "#1F263C", "#36574E",
            "#52CE79", "#ADAAA9", "#8A9F45", "#6542D2", "#00FB8C", "#5D697B", "#CCD27F", "#94A5A1",
            "#790229", "#E383E6", "#7EA4C1", "#4E4452", "#4B2C00", "#620B70", "#314C1E", "#874AA6",
            "#E30091", "#66460A", "#EB9A8B", "#EAC3A3", "#98EAB3", "#AB9180", "#B8552F", "#1A2B2F",
            "#94DDC5", "#9D8C76", "#9C8333", "#94A9C9", "#392935", "#8C675E", "#CCE93A", "#917100",
            "#01400B", "#449896", "#1CA370", "#E08DA7", "#8B4A4E", "#667776", "#4692AD", "#67BDA8",
            "#69255C", "#D3BFFF", "#4A5132", "#7E9285", "#77733C", "#E7A0CC", "#51A288", "#2C656A",
            "#4D5C5E", "#C9403A", "#DDD7F3", "#005844", "#B4A200", "#488F69", "#858182", "#D4E9B9",
            "#3D7397", "#CAE8CE", "#D60034", "#AA6746", "#9E5585", "#BA6200")
######################## SET UP LINES TO RUN ########################

set_cores <- function(x){
  CORES_TO_USE <- x
  registerDoMC(cores=CORES_TO_USE)
}

registerDoMC(cores=CORES_TO_USE)

#####################################################################
########## Main Pipeline Code to run RF classifcation ###############
#####################################################################

### function takes in a feature table containing all features that
### we want to use to for classifcation and a second vector with the
### expected classes, it will also take in a metric to train and save
### results for as well as a sampling procedure
### it is expected that features have already been filtered to remove
### uncommonly found features also note that by default this function
### splits data 80-20 and runs cross validation
rf_classification_pipeline <- function(feature_table, classes, metric, 
                                      ntree, nmtry, sampling,
                                      nfolds, ncrossrepeats, pro, SEED){
  #create vectors to save cv and test metric values for each data split 
  message("Make sure that the disease class is level 1 (called Case), and the control (called Control) is level 2!!!
          if this is not the case results from this function will be wrong for PR")
  
  
  ### create data split
  message("splitting data")
  
  ### Head function will need to pass in SEED values to run 100 data partitions or it won't be reproducible
  set.seed(SEED)
  train_index <- createDataPartition(classes, p=pro, list=FALSE)
  message(train_index)
  ## make training feature and class tables
  train_classes <- classes[train_index]
  train_features <- feature_table[train_index,]
  ## make test feature and class tables
  test_classes <- classes[-train_index]
  test_features <- feature_table[-train_index,]
  message("finished splitting data")
  ### do we want to do cross fold validation or montecarlo??
  ### for now we will default to cv
  
  folds <- nfolds
  repeats <- ncrossrepeats
  mtry_num <- nmtry
  ## change back to 100 in future
  set.seed(SEED)
  cvIndex <- caret::createMultiFolds(factor(train_classes), folds, times=repeats)
  set.seed(SEED)
  seeds_len <- (repeats*folds) + 1
  seeds <- vector(mode="list", length=(seeds_len))
  for(i in 1:(seeds_len-1)){
    seeds[[i]] <- sample.int(n=1000, mtry_num)
  }
  seeds[seeds_len] <- sample.int(1000, 1)
  if(metric=="PR"){
    message("Metric to test on is AUPRC")
    cv <- trainControl(method="repeatedcv",
                       number=folds,
                       index=cvIndex,
                       returnResamp = "final",
                       summaryFunction=prSummary,
                       classProbs=TRUE,
                       savePredictions = TRUE,
                       seeds=seeds)
    metric="AUC"
    }else if(metric=="ROC"){
      cv <- trainControl(method="repeatedcv",
                         number=folds,
                         index=cvIndex,
                         returnResamp = "final",
                         classProbs=TRUE,
                         summaryFunction=twoClassSummary,
                         savePredictions = TRUE,
                         seeds=seeds)
      metric <- "ROC"
    }
  ### set sampling method
  cv$sampling <- sampling
  ### set up grind to search for best mtry
  
  ### set up mtry values to test
  n_features <- ncol(feature_table)
  message(n_features)
  #set up 6 different mtry values to test
  mtry <- round(seq(1, n_features/3, length=nmtry))
  mtry <- mtry[mtry <= n_features]
  message(mtry)
  grid <- expand.grid(mtry = mtry)
  ## train the model
  message("Training model")
  set.seed(SEED)
  trained_model <- train(train_features, train_classes,
                         method="rf",
                         trControl=cv,
                         metric=metric,
                         tuneGrid=grid,
                         ntree=ntree,
                         importance=TRUE)
  
  message("Finished training model")
  # Mean AUC value over repeates of the best hyperparameter during training
  if(metric=="ROC"){
    cv_auc <- getTrainPerf(trained_model)$TrainROC
  }else if(metric=="AUC"){
    cv_auc <- getTrainPerf(trained_model)$TrainAUC
  }
  
  ### get important features for later validation...
  important_features <- trained_model$finalModel$importance
  
  ## predict on the test set and get predicted probabilities
  rpartProbs <- predict(trained_model, test_features, type="prob")
  if(metric=="ROC"){
    ### critcal that factor levels are correct for this to calculate.... 
    test_roc <- pROC::roc(ifelse(test_classes=="Case", 1, 0), rpartProbs[[2]])
    test_auc <- test_roc$auc
  }else if(metric=="AUC"){
    ######### this doesn't compute right needs to be fixed!!!!
    #get probs for postive class
    matriz <- cbind(test_classes, predict(trained_model, test_features, type="prob"), predict(trained_model, test_features))
    names(matriz) <- c("obs", levels(test_classes), "pred")
    pr_test_stats <- prSummary(matriz, levels(test_classes))
    test_auc <- pr_test_stats[1]
  }else{
    message("metric not used by this pipeline")
    return(NULL)
  }

  ### return in a list the cv_auc, the test_auc, the trained_model results and the important features...
  results <- list(cv_auc, test_auc, trained_model$results, important_features, trained_model)
  return(results)
}
#####################################################################
############# Main Pipeline to run RF regression ####################
rf_regression_pipeline <- function(feature_table, actual, SEED, sampling ){
  
  message("Running random forest regression pipeline")
  message("Make sure that the actual variable is a numeric class")
  
  message("Splitting data into test and train")
  set.seed(SEED)
  #setting indexs for training and test data
  train_index <- createDataPartition(actual, p=.8, list=FALSE)
  train_features <- feature_table[train_index,]
  test_features <- feature_table[-train_index,]
  train_actual <- actual[train_index]
  test_actual <- actual[-train_index]
  
  ### okay data is split
  folds <- 5
  repeats <- 10
  mtry_num <- 7
  
  set.seed(SEED)
  cvIndex <- createMultiFolds(train_actual, folds, times=repeats)
  set.seed(SEED)
  seeds_len <- (repeats*folds) + 1
  seeds <- vector(mode="list", length=(seeds_len))
  for(i in 1:(seeds_len-1)){
    seeds[[i]] <- sample.int(n=1000, mtry_num)
  }
  seeds[seeds_len] <- sample.int(1000, 1)
  cv <- trainControl(method="repeatedcv",
                     number=folds,
                     index=cvIndex,
                     returnResamp = "final",
                     classProbs = FALSE,
                     indexFinal = NULL,
                     savePredictions = TRUE,
                     seeds = seeds)
  cv$sampling <- sampling
  ### set up grind to search for best mtry
  
  ### set up mtry values to test
  n_features <- ncol(feature_table)
  
  #set up 6 different mtry values to test
  mtry <- floor(seq(1, n_features/3, length=6))
  mtry <- mtry[mtry <= n_features]
  mtry <- c(mtry, 2)
  grid <- expand.grid(mtry = mtry)
  message(mtry)
  message("Training model")
  set.seed(SEED)
  trained_model <- train(train_features, train_actual,
                         method="rf",
                         trControl=cv,
                         metric="RMSE",
                         tuneGrid=grid,
                         ntree=2001,
                         importance=TRUE)
  #get RMSE value over repeats for the best mtry parameter
  cv_best <- getTrainPerf(trained_model)
  cv_results <- trained_model$results
  
  #get RMSE value for best model on test daata
  predictions <- predict(trained_model, test_features, type="raw")
  test_results <- postResample(test_actual, predictions)
  
  #get important features
  important_features <- trained_model$finalModel$importance
  
  results <- list(cv_best, cv_results, test_results, important_features,
                  trained_model)
  return(results)
}


################## Function to run RF pipelines #####################

############################## RF CLASSIFCATION RUN #################
#### main function use to run RF with will produce a dataframe and save
### that dataframe to the path that is specified
get_rf_results <- function(feature_table, classes, metric="ROC", sampling=NULL, 
                           repeats=10, path, nmtry=6, ntree=1001, 
                           nfolds=3, ncrossrepeats=10, pro=0.8, list_of_seeds, help=FALSE,
                           version=FALSE, values=FALSE){
  if(help==TRUE){
    message("Example of function usage: get_rf_results(feature_table = input_features, classes = classes, metric='ROC')")
    message("Explanation of each parameter:")
    message("feature_table: The input feature table that you want to use to predict the input to classes. The format of this
            table should be that rows equal samples and columns equal features (such as abundance of each microbe).")
    message("classes: The classes of each sample that is found in the rows of the feature table. Make sure that the classes
            order matches the sample order in feature_table.")
    message("metric: This is the metric that you want to evaluate your model with. The default for this is 'ROC' which will
            evaluate model performance based on  AUC of a ROC curve. Other available options are 'PR' which evaluates 
            based on the AUC of a precision recall curve, this is useful for unbalanced class designs.")
    message("sampling: This is used to set the sampling strategy during model training. The default is NULL meaning no
            special sampling will be done and should be used in balance class designs. Options include up-sampling, 
            down-sampling, and smote-sampling.")
    message("repeats: the number of training and test datasplits you want to be done. Default: 10")
    message("path: this is the path that you want to save the results of this pipeline to.")
    message("nmtry: is the number of mtry parameters you want to test during model optimisation. Default: 3")
    message("ntree: the number of trees you want you model to use. Default: 1001")
    message("nfolds: the number of even folds to split your data into during cross validation. Default: 3")
    message("ncrossrepeats: the number of times you want to repeat cross validation for each data split. Default: 10")
    message("pro: The proportion of data that is used as training data. The other proportion will be used as a test
            /validation dataset. Default: 0.8")
    message("list_of_seeds: A list of random seeds that must be equal to or longer than the number input into repeats.")
    message("help: Set to TRUE in order to see these messages!")
    message("version: set to TRUE to print out the version of this pipeline")
    message("values: set to TRUE to print out the expected results from this function")
    return(NULL)
  }else if(values==TRUE){
    message("This function returns an object with the follow characteristics:")
    message("Object[[1]] contains all the median cross validation AUCS from each data split using the best mtry value")
    message("Object[[2]] contains all the test AUC values from each data split")
    message("Object[[3]] contains all the tested mtry values and the median ROC for each from each data split")
    message("Object[[4]] contains the list of important features from the best model selected from each data split")
    message("Object[[5]] contains each caret random forest model from each data split")
    message("This function will also write a csv with cv AUCS and test AUCS, to the given path as well as an RDS file that
            contains the resulting object from this function")
    return(NULL)
  }
  if(version==TRUE){
    message("Version 0.1.1")
  }
  start_time <- Sys.time()
  cv_aucs <- c()
  test_aucs <- c()
  results_total <- list()
  important_features <- list()
  models <- list()
  #run rf for the number of repeats given
  for(i in 1:repeats){
    message("On Round ",i)
    int_res <- rf_classification_pipeline(feature_table = feature_table, classes = classes, 
                           metric = metric, sampling = sampling, 
                           SEED = list_of_seeds[i],
                           nmtry= nmtry,
                           ntree=ntree,
                           nfolds=nfolds,
                           ncrossrepeats=ncrossrepeats,
                           pro = pro)
    cv_aucs <- c(cv_aucs, int_res[[1]])
    test_aucs <- c(test_aucs, int_res[[2]])
    list_names <- paste0("t",i)
    results_total[[list_names]] <- int_res[[3]]
    important_features[[list_names]] <- int_res[[4]]
    colnames(important_features[[list_names]]) <- gsub("^",list_names,colnames(important_features[[list_names]]))
    models[[list_names]] <- int_res[[5]]
  }
  #take data save make one master list to save it all as an rds... (given the path)
  #take the test_auc and cv_auc and write it out into a csv
  #take the results list rbind thw whole thing and write it out as a csv
  #take the important features list rbind the whole thing and write it out as a csv
  #don't do anything with models just return it in the master list and save the rds as above
  
  #master ret list
  ret_list <- list(cv_aucs, test_aucs, results_total, important_features,
                   models)
  
  #auc dataframe
  auc_data <- data.frame(cv_auc = cv_aucs, test_auc = test_aucs)
  write.csv(auc_data, file=paste0(path,"_aucs.csv"))
  
  #take results list and rbind it
  cv_results <- bind_rows(results_total)
  write.csv(cv_results, file=paste0(path,"_cv_results.csv"))
  
  ret_features <- do.call(cbind, important_features)
  write.csv(ret_features, file=paste0(path,"_imprt_feats.csv"))

  #save rds
  saveRDS(ret_list, file=paste0(path,"_masterlist.rds"))
  endtime <- Sys.time()
  total_time <- endtime - start_time
  message(total_time)
  return(ret_list)
}

#################### Random RF Testing ##############################
get_random_rf_results <- function(feature_table, list_of_scrambles, metric="ROC", sampling=NULL, 
                                  repeats=10, path, nmtry=6, ntree=1001, 
                                  nfolds=3, ncrossrepeats=10, pro=0.8, list_of_seeds, help=FALSE,
                                  version=FALSE, values=FALSE){
  message("Number of repeats must be equal to number of scramble lists")
  
  if(help==TRUE){
    message("Example of function usage: get_rf_results(feature_table = input_features, classes = classes, metric='ROC')")
    message("Explanation of each parameter:")
    message("feature_table: The input feature table that you want to use to predict the input to classes. The format of this
            table should be that rows equal samples and columns equal features (such as abundance of each microbe).")
    message("list_of_scrambles: Randomly scrambled classifcations for each sample. Not that a new randomization is used at
            each data split and so the length of this list must be equal to the number of repeats")
    message("metric: This is the metric that you want to evaluate your model with. The default for this is 'ROC' which will
            evaluate model performance based on  AUC of a ROC curve. Other available options are 'PR' which evaluates 
            based on the AUC of a precision recall curve, this is useful for unbalanced class designs.")
    message("sampling: This is used to set the sampling strategy during model training. The default is NULL meaning no
            special sampling will be done and should be used in balance class designs. Options include up-sampling, 
            down-sampling, and smote-sampling.")
    message("repeats: the number of training and test datasplits you want to be done. Default: 10")
    message("path: this is the path that you want to save the results of this pipeline to.")
    message("nmtry: is the number of mtry parameters you want to test during model optimisation. Default: 3")
    message("ntree: the number of trees you want you model to use. Default: 1001")
    message("nfolds: the number of even folds to split your data into during cross validation. Default: 3")
    message("ncrossrepeats: the number of times you want to repeat cross validation for each data split. Default: 10")
    message("pro: The proportion of data that is used as training data. The other proportion will be used as a test
            /validation dataset. Default: 0.8")
    message("list_of_seeds: A list of random seeds that must be equal to or longer than the number input into repeats.")
    message("help: Set to TRUE in order to see these messages!")
    message("version: set to TRUE to print out the version of this pipeline")
    message("values: set to TRUE to print out the expected results from this function")
    return(NULL)
  }else if(values==TRUE){
    message("This function returns an object with the follow characteristics:")
    message("Object[[1]] contains all the median cross validation AUCS from each data split using the best mtry value")
    message("Object[[2]] contains all the test AUC values from each data split")
    message("Object[[3]] contains all the tested mtry values and the median ROC for each from each data split")
    message("Object[[4]] contains the list of important features from the best model selected from each data split")
    message("Object[[5]] contains each caret random forest model from each data split")
    message("This function will also write a csv with cv AUCS and test AUCS, to the given path as well as an RDS file that
            contains the resulting object from this function")
    return(NULL)
  }
  if(version==TRUE){
    message("Version 0.1")
  }
                                    
                                                          
  
  start_time <- Sys.time()
  cv_aucs <- c()
  test_aucs <- c()
  results_total <- list()
  important_features <- list()
  models <- list()
  #run rf for the number of repeats given
  for(i in 1:repeats){
    message("On Round ",i)
    int_res <- rf_classification_pipeline(feature_table = feature_table, 
                                          classes = list_of_scrambles[[i]], 
                                          metric = metric, sampling = sampling, 
                                          SEED = list_of_seeds[i],
                                          nmtry= nmtry,
                                          ntree=ntree,
                                          nfolds=nfolds,
                                          ncrossrepeats=ncrossrepeats,
                                          pro = pro)
    cv_aucs <- c(cv_aucs, int_res[[1]])
    test_aucs <- c(test_aucs, int_res[[2]])
    list_names <- paste0("t",i)
    results_total[[list_names]] <- int_res[[3]]
    important_features[[list_names]] <- int_res[[4]]
    colnames(important_features[[list_names]]) <- gsub("^",list_names,colnames(important_features[[list_names]]))
    models[[list_names]] <- int_res[[5]]
  }
  #take data save make one master list to save it all as an rds... (given the path)
  #take the test_auc and cv_auc and write it out into a csv
  #take the results list rbind thw whole thing and write it out as a csv
  #take the important features list rbind the whole thing and write it out as a csv
  #don't do anything with models just return it in the master list and save the rds as above
  
  #master ret list
  ret_list <- list(cv_aucs, test_aucs, results_total, important_features,
                   models)
  
  #auc dataframe
  auc_data <- data.frame(cv_auc = cv_aucs, test_auc = test_aucs)
  write.csv(auc_data, file=paste0(path,"_aucs.csv"))
  
  #take results list and rbind it
  cv_results <- bind_rows(results_total)
  write.csv(cv_results, file=paste0(path,"_cv_results.csv"))
  
  ret_features <- do.call(cbind, important_features)
  write.csv(ret_features, file=paste0(path,"_imprt_feats.csv"))
  
  #save rds
  saveRDS(ret_list, file=paste0(path,"_masterlist.rds"))
  endtime <- Sys.time()
  total_time <- endtime - start_time
  message(total_time)
  return(ret_list)
}
##################### RF Regression Run #############################
get_rf_regres_results <- function(feature_table, actual, sampling,
                                  repeats, path, list_of_seeds){
  start_time <- Sys.time()
  cv_best <- list()
  cv_results <-list()
  test_results <- list()
  important_features <- list()
  models <- list()
  
  for(i in 1:repeats){
    message("On Round ", i)
    int_res <- rf_regression_pipeline(feature_table, actual,
                                      sampling=sampling, 
                                      SEED=list_of_seeds[i])
    list_name <- paste0("t",i)
    cv_best[[list_name]] <- int_res[[1]]
    cv_results[[list_name]] <- int_res[[2]]
    test_results[[list_name]] <- int_res[[3]]
    important_features[[list_name]] <- int_res[[4]]
    colnames(important_features[[list_name]]) <- 
      gsub("^",list_name,colnames(important_features[[list_name]]))
    models[[list_name]] <- int_res[[5]]
  }
  

  cv_best_data <- bind_rows(cv_best)
  write.csv(cv_best_data, file=paste0(path, "_cv_best.csv"))
  
  cv_results_data <- bind_rows(cv_results)
  write.csv(cv_results_data, file=paste0(path, "_cv_results.csv"))
  
  test_results <- bind_rows(test_results)
  write.csv(test_results, file=paste0(path,"_test_results.csv"))
  
  impt_feats <- do.call(cbind, important_features)
  write.csv(impt_feats, file=paste0(path,"_impt_feats.csv"))

  
  
  end_time <- Sys.time()
  totaltime <- end_time - start_time
  ret_list <- list(cv_best, cv_results, test_results,
                   important_features, models, totaltime)
  saveRDS(ret_list, file = paste0(path,"_master_list.rds"))  
  return(ret_list)
}


#####################################################################
###################### Utility Functions ############################
#####################################################################
Calc_mean_accuray_decrease <- function(impt_feat_list, help=FALSE){
  
  if(help==TRUE){
    message("This function takes in a list of dataframes that contain the important
            feature stats for each random forest training run")
    message("It will return a dataframe containing the mean, sd, max, and min decrease
            in accuracy for each feature over all of the random forest training rounds")
  }
  
  pass=TRUE
  #make sure colnames for each impt_feat_table matches up
  for(j in 1:length(impt_feat_list)){
    if(!all.equal(rownames(impt_feat_list[[1]]), rownames(impt_feat_list[[i]]))){
      pass=FALSE
    }
  }
  
  if(pass==FALSE){
    message("rownames don't match")
    return(NULL)
  }
  
  ret_data_list <- list()
  for(i in 1:length(impt_feat_list)){
    prefix <- paste0("t",i)
    message(prefix)
    acc <- impt_feat_list[[i]][,3]
    ret_data_list[[prefix]] <- acc
  }
  
  
  ret_data_frame <- do.call(cbind, ret_data_list)
  test <- as.data.frame(cbind(ret_data_frame, t(apply(ret_data_frame, 1, Stats))))
  ret_final_frame <- test[order(test$Mean, decreasing=TRUE), ]
  return(ret_final_frame)
}

#####################################################################
##################### Stats from Rows ###############################
Stats <- function(x){
  Mean <- mean(x, na.rm=TRUE)
  SD <- sd(x, na.rm=TRUE)
  Min <- min(x, na.rm=TRUE)
  Max <- max(x, na.rm=TRUE)
  return(c(Mean=Mean, SD=SD, Min=Min, Max=Max))
}



#####################################################################
### Generate ROC curves

generate_ROC_curve <- function(RF_models, dataset, labels, title, help=F){
  
  if(help){
    
    message("RF_model: a list of RF models generated by caret. If using the pipeline from this script they will be contained in the fifth index of the returned list")
    message("dataset: the input features that were using to create the RF models")
    message("labels: a dataframe with the column classes that defines the class labels for each input sample")
    message("title: the title at the top of the generated plot")
    message("help: set this true to see this message!")
  }
  
  AUC_data <- vector()
  ROC_data <- list()
  ROC_curve_data <- list()
  
  ROC_sens <- list()
  
  ROC_specs <- list()
  
  for(i in 1:length(RF_models)){
    
    prefix <- paste0("t",i)
    #grab the samples that the model was trained on
    training_samples <- rownames(RF_models[[i]]$trainingData)
    
    
    #get the samples to predict data from
    prediction_set <- dataset[!rownames(dataset) %in% training_samples,]
    prediction_set_labels <- labels[!rownames(labels) %in% training_samples,,drop=F]
    message(prediction_set_labels$classes)
    message("getting hold-out data from each model training session for test validation")
    
    #make predictions on this dataset using the final model
    
    
    predicts <- predict(RF_models[[i]], prediction_set, type="prob")
    message("making predictions")
    message(length(prediction_set_labels))
    roc_data <- pROC::roc(prediction_set_labels$classes, predicts[,1], levels=c("Control", "Case"))
    
    sens <- roc_data$sensitivities
    specs <- roc_data$specificities
    
    indexs_to_keep <- floor(seq(1, length(sens), length=40))
    ROC_sens[[prefix]] <- sens[indexs_to_keep]
    ROC_specs[[prefix]] <- specs[indexs_to_keep]
    AUC_data <- c(AUC_data, roc_data$auc)
    if(i==1){
      plot(roc_data, xlim=c(1,0), ylim=c(0,1), col="red")
    }else{
      plot(roc_data,add=T, xlim=c(1,0), ylim=c(0,1))
    }
  }
  
  ROC_curve_data[[1]] <- do.call(cbind, ROC_sens)
  ROC_curve_data[[2]] <- do.call(cbind, ROC_specs)
  
  #turn data into long form
  
  SENS_melt <- melt(ROC_curve_data[[1]])
  SPEC_melt <- melt(ROC_curve_data[[2]])
  
  if(all.equal(SENS_melt$Var1, SPEC_melt$Var1) & all.equal(SENS_melt$Var2, SPEC_melt$Var2)){
    SENS_melt$Value2 <- SPEC_melt$value
  }else(
    message("values dont match")
  )
  #alright so we have the sens and spec for each thing so lets generate curves from each column
  
  #calc mean values this part i'm not 100% sure how to do....
  
  
  SENs_values <- as.data.frame(cbind(ROC_curve_data[[1]], t(apply(ROC_curve_data[[1]], 1, Stats))))
  SPEC_values <- as.data.frame(cbind(ROC_curve_data[[2]], t(apply(ROC_curve_data[[2]], 1, Stats))))
  
  mean_values <- data.frame(Sens= SENs_values$Mean,
                            Specs= SPEC_values$Mean)
  
  AUC_label <- paste("Mean AUC", round(mean(AUC_data), digits=5))
  plot <- ggplot() + geom_path(data=SENS_melt, mapping=aes(x=Value2, y=value, group=Var2), alpha=0.05) + xlab("Specificity") +
    ylab("Sensitivity") + scale_x_reverse() + scale_color_manual(values = COLORS) + theme(legend.position = "none") + 
    geom_abline(intercept=1, slope=1) + geom_line(data=mean_values, aes(x=Specs, y=Sens), color="Red") + ggtitle(title) +
    annotate(geom="text", y=0, x=0.2, label=AUC_label, color="red")
  
  return(plot)
}




