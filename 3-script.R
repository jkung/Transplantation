##################### 
# Load in the data
# Select OPOs and ABOs to run for
#####################

## First make sure to set the working directory properly:
setwd(".")

load("kidney.Rdata")

source("0-functions.R")
source("1-timings.R")

# list of the top ten OPOs by population
OPOs <- c("CAOP-OP1", "TXGC-OP1", "TXSB-OP1", "NYRT-OP1", "PADV-OP1","ILIP-OP1", "CADN-OP1", "MAOB-OP1", "MIOP-OP1", "NCNC-OP1")
ABOs <- c("O", "A", "B", "AB")

# toggle between 1 and 38 to run for all time configurations, or just the one with most data
# 32 = 3m with most data; 34 = 6m with most data; 38 = 12m with most data
times = 1:38

desiredBucket = 2 # the kidney quality that we are willing to accept (1 is <= 0.2; 2 is <= 0.4; etc.)
# this can be toggled for other kidney quality thresholds

ID = sprintf("%s-Q%d.Rdata", paste(OPOs, collapse = "-"), desiredBucket)

##################### 
# Restrict to OPOs above, and to non-bypassed
#####################
PTRdata <- PTRdata %>% tbl_df
PTRcopy <- PTRdata
PTRdata <- PTRdata %>% filter(LISTING_CTR_DSA %in% OPOs,
                              OPO_CTR %in% OPOs,
                              BYPASSED == 0,
                              !is.na(KDPI_CALC),
                              ORGAN_SEQUENCE_NUM <= 500) %>%
  mutate(MATCH_DATE = as.POSIXct(as.character(MATCH_DATE), format = "%Y-%m-%d"),
         INIT_DATE = as.POSIXct(as.character(INIT_DATE), format = "%Y-%m-%d"),
         daysSince = difftime(MATCH_DATE, INIT_DATE, units = "days"),
         smallDR = pmin(DR1, DR2),
         largeDR = pmax(DR1, DR2),
         DR1 = smallDR, 
         DR2 = largeDR,
         kidneyquality = ((KDPI_CALC-.01) %/% 0.2 + 1)) %>% #0-.2 inclusive = 1, 0.21-0.4 inclusive = 2, etc.
  select(-smallDR, -largeDR)

sprintf("%d rows in PTR originally", nrow(PTRcopy))
sprintf("%d rows removed after restricting to OPOs, no bypasses, non-NA KDPI", nrow(PTRcopy) - nrow(PTRdata))
sprintf("%d rows remain in PTR", nrow(PTRdata))

KidneyWLcopy <- KidneyWaitlist

sprintf("%d rows originally in KidneyWaitlist", nrow(KidneyWLcopy))

##################### 
# Remove directed donations
#####################
PTRdata <- PTRdata %>%
  mutate(directedText = as.numeric(grepl("DIRECT", REFUSAL_OSTXT, ignore.case = TRUE)))
kidneysToRemove <- PTRdata %>% group_by(DONOR_ID) %>%
  summarize(allRefusalReasons = paste(unique(REFUSAL_ID), collapse = "*", sep = ""),
            directed = grepl("851", allRefusalReasons),
            dirSum = sum(directedText),
            toRemove = (directed | dirSum)) 
removeDonorID <- kidneysToRemove %>% filter(toRemove) %>% .$DONOR_ID %>% unique
PTRdata <- PTRdata %>% filter(!(DONOR_ID %in% removeDonorID))

sprintf("%d kidneys removed due to directed donation", sum(kidneysToRemove$toRemove == TRUE))

temp <- getPatientsKidneys()
patients <- temp[[1]]
donors <- temp[[2]]

firstAccept <- PTRdata %>% group_by(DONOR_ID) %>%
  summarize(firstAccept = ORGAN_SEQUENCE_NUM[ACCEPTED == 1][1])
donors <- donors %>% inner_join(firstAccept)

sprintf("%d patients in PTR", nrow(patients))
sprintf("%d donors in PTR", nrow(donors))
rm(temp)

##################### 
# Clean up the KidneyWaitlist data 
#####################
KidneyWaitlist <- addNextDate()
KidneyWaitlist <- filterDates("2000-05-01", "2013-06-30") 
KidneyWaitlist <- KidneyWaitlist %>% 
  mutate(CHG_DATE = as.POSIXct(as.character(CHG_DATE), format = "%Y-%m-%d"),
         nextDate = as.POSIXct(as.character(nextDate), format = "%Y-%m-%d"),
         totalDays = as.numeric(difftime(nextDate, CHG_DATE, units = "days")))

##################### 
# Create data frame of all patients, not just ones in PTR
# Also create the expected DR mismatch table
#####################
load("KidpanSubset.RData")
load("patientlocs.Rdata")

patientlocs <- patientlocs %>% 
  mutate(opo = substring(LISTING_CTR_DSA, 1, 8)) %>% select(-LISTING_CTR_DSA)

KidpanSubset <- KidpanSubset %>%
  mutate(ABO = CAND_ABO_AT_REGISTRATION) %>%
  select(-CAND_ABO_AT_REGISTRATION, -LISTING_CTR_DSA)

intermed <- KidpanSubset %>% inner_join(patientlocs) %>% 
  filter(opo %in% OPOs,
         ABO %in% ABOs)

noOfferPatients <- intermed %>% filter(!(WL_ID_CODE %in% patients$WL_ID_CODE)) %>%
  mutate(numOffers = 0,
         numAccepted = 0,
         initdate = as.POSIXct(as.character(INIT_DATE, format = "%Y-%m-%d")),
         opo = opo,
         age = INIT_AGE,
         dr1 = pmin(DR1, DR2),
         dr2 = pmax(DR1, DR2),
         abo = ABO,
         lastmatch = as.POSIXct(as.character(END_DATE, format = "%Y-%m-%d"))) %>%
  select(WL_ID_CODE, initdate, numAccepted, opo, age, dr1, dr2, numOffers, abo, lastmatch) %>% tbl_df %>%
  filter(lastmatch >= as.POSIXct("2007-07-01", format = "%Y-%m-%d"),
         initdate >= as.POSIXct("2000-01-01", format = "%Y-%m-%d"))

allPatientsPre <- patients %>% rbind(noOfferPatients) 
  
sprintf("%d patients without any offers added to patient list", nrow(noOfferPatients))

##################### 
# Set the parameters here (prediction range, amount of train/valid/test data, OPO, blood group)
#####################
#varies for each time frame we want to predict on

timingIndex = NA
predTime = NA
OPO = NA
ABO = NA
nTrainVal = NA
nTrainValTrue = NA
modelIndex = NA
nTest = NA
nTestTrue = NA
log.valid.auc = NA
log.test.auc = NA
cart.test.auc = NA
rf.valid.auc = NA
rf.test.auc = NA

DRmatching <- read.table("kidpan_data.dat",sep=',', quote = "\"",header=TRUE) %>% 
  select(WL_ID_CODE, DR1, DR2, DONOR_ID, DDR1, DDR2, DRMIS)

df = data.frame(timingIndex, predTime, OPO, ABO,
                nTrainVal, nTrainValTrue,
                modelIndex,
                nTest, nTestTrue,
                log.valid.auc, log.test.auc,
                cart.test.auc, 
                rf.valid.auc, rf.test.auc)

df.var.Imp <- data.frame() %>% 
  mutate(opo = NA,
         abo = NA,
         q12intensity = NA,
         medIntensity = NA,
         lowIntensity = NA,
         avgDRMISprev = NA,
         avgOSQprev = NA,
         avgKDPIprev = NA,
         avgMIS = NA,
         cpra = NA,
         cpraFactor = NA,
         years = NA,
         pediatric = NA,
         month = NA,
         timingInd = NA)

for (timeIndex in times){
  print(timeIndex)
  i <- timeIndex
  
  predRangeAmount = timings[i, "predRange"]
  trainDataAmount = timings[i, "trainData"]
  validDataAmount = timings[i, "validData"]
  testDataAmount = timings[i, "testData"]
  
  TRAINING_START_DATE = as.POSIXct("2007-05-01", format = "%Y-%m-%d")
  TRAINING_END_DATE = TRAINING_START_DATE + days(trainDataAmount)
  
  VALIDATION_START_DATE = TRAINING_END_DATE + days(predRangeAmount) + days(1)
  VALIDATION_END_DATE = VALIDATION_START_DATE + days(validDataAmount)
  
  TESTING_START_DATE = VALIDATION_END_DATE + days(predRangeAmount) + days(1)
  TESTING_END_DATE = TESTING_START_DATE + days(testDataAmount)
  
  pivotDates = c(TRAINING_START_DATE, TRAINING_END_DATE, VALIDATION_START_DATE, VALIDATION_END_DATE, TESTING_START_DATE, TESTING_END_DATE)
  
  source("2-DRmatching.R")
  ##################### 
  # Create each of the training/validation/testing sets
  #####################
  
  firstPoints <- pivotDates[c(1, 3, 5)]
  lastPoints <- pivotDates[c(2, 4, 6)]
  
  firstDates <- KidneyWLcopy %>% 
    filter(WL_ID_CODE %in% allPatients$WL_ID_CODE) %>%
    group_by(WL_ID_CODE) %>% 
    summarize(firstDate = first(CHG_DATE)) %>%
    mutate(firstDate = as.POSIXct(as.character(firstDate), format = "%m/%d/%Y")) %>%
    select(WL_ID_CODE, firstDate) %>%
    ungroup()
  
  lastDates <- KidneyWaitlist %>%
    select(WL_ID_CODE, CHG_DATE, nextDate, UNOS_CAND_STAT_CD, totalDays) %>%
    group_by(WL_ID_CODE) %>%
    summarize(lastDate = as.POSIXct(as.character(last(nextDate)), format = "%Y-%m-%d"))
  
  outSet <- vector("list", length = 3)
  
  for (index in seq(1, 3))
  {
    dates <- lastDates %>% inner_join(firstDates)
    
    firstBoundary = dates$firstDate
    firstBoundary[firstBoundary <= firstPoints[index]] = firstPoints[index]
    dates$firstDate = firstBoundary
    
    lastBoundary = dates$lastDate
    lastBoundary[lastBoundary >= lastPoints[index]] = lastPoints[index]
    dates$lastDate = lastBoundary
    
    dates <- dates %>% mutate(num3months = as.numeric(
      difftime(lastDate, firstDate, units = "weeks") / 25) %>% round %>% '*'(4)) %>%
      filter(num3months >= 0)
    # this corresponds to a sampling rate of 8 observations per patient per year
    
    set.seed(81816)
    
    outputSet <- rep(dates$WL_ID_CODE, dates$num3months)
    outputSet = data.frame(outputSet)
    names(outputSet) = "WL_ID_CODE"
    outputSet$predictDate = sample(seq(firstPoints[index], lastPoints[index], by = "day"),
                                   nrow(outputSet), replace = TRUE) %>% 
      substring(1, 10) %>% 
      as.character %>% 
      as.POSIXct(format = "%Y-%m-%d")
    
    outSet[[index]] <- outputSet %>% inner_join(dates) %>%
      filter(predictDate >= firstDate,
             predictDate <= lastDate)
  }
  
  outputSet = rbind(outSet[[1]], outSet[[2]], outSet[[3]])
  
  nrow(allPatients)
  nrow(patients)
  nrow(noOfferPatients)
  
  merged <- outputSet %>% select(-lastDate, -firstDate, -num3months) %>%
    mutate(uniqueID = row_number()) %>% tbl_df
  merged2 <- merged
  merged3 <- merged # merged4 in old is merged3 here
  merged4 <- merged # temp in old is merged 4 here
  original <- merged
  
  ##################### 
  # Compute some of the independent variables
  #####################
  
  merged <- merged %>%
    inner_join(allPatients %>%
                 select(initdate, opo, age, abo, avgMIS, WL_ID_CODE)) %>%
    inner_join(PTRdata %>%
                 select(WL_ID_CODE, ORGAN_SEQUENCE_NUM,
                        MATCH_DATE, DRMIS, CDC_RISK_HIV_DON,
                        HEP_C_ANTI_DON, KIDNEY_CATEGORY, KDPI_CALC, ACCEPTED, 
                        kidneyquality)) %>%
    mutate(isPrev = MATCH_DATE < predictDate,
           isUnder50 = ORGAN_SEQUENCE_NUM <= 50,
           toSum = 1,
           predictEnd = predictDate + days(predRangeAmount)) %>%
    group_by(uniqueID) %>%
    summarize(numSeen = sum(isPrev),
              numq12seen = sum(toSum[isPrev & kidneyquality <= 2]),
              numMediumSeen = sum(toSum[isPrev & (KDPI_CALC > .40 & KDPI_CALC <= .70)]),
              numLowSeen = sum(toSum[isPrev & KDPI_CALC > .70]),
              avgDRMISprev = mean(DRMIS[isPrev]),
              avgOSQprev = median(ORGAN_SEQUENCE_NUM[isPrev]),
              avgKDPIprev = mean(KDPI_CALC[isPrev], na.rm = TRUE),
              acceptsNext = sum(ACCEPTED[MATCH_DATE < predictEnd]),
              daysUntilDesired = min(as.numeric(
                difftime(MATCH_DATE[kidneyquality <= desiredBucket &
                                      !isPrev &
                                      isUnder50],
                         predictDate,
                         units = "days"))))
  
  mergedNonOffers <- original %>% filter(WL_ID_CODE %in% noOfferPatients$WL_ID_CODE) %>%
    mutate(numSeen = 0,
           numq12seen = 0,
           numMediumSeen = 0,
           numLowSeen = 0,
           avgDRMISprev = NA,
           avgOSQprev = NA,
           avgKDPIprev = NA,
           acceptsNext = 0,
           daysUntilDesired = Inf) %>%
    select(uniqueID,
           numSeen,
           numq12seen,
           numMediumSeen,
           numLowSeen,
           avgDRMISprev,
           avgOSQprev,
           avgKDPIprev,
           acceptsNext,
           daysUntilDesired)
  
  merged <- rbind(merged, mergedNonOffers)
  
  ##################### 
  # Compute current cpra
  #####################
  
  cpraWL <- KidneyWLcopy %>% select(WL_ID_CODE, CHG_TY, CHG_DATE, CHG_TIME,
                                    CPRA) %>% tbl_df %>%
    mutate(CHG_DATE = as.POSIXct(CHG_DATE, format = "%m/%d/%Y"),
           CPRA = round(CPRA, 2)) %>%
    group_by(WL_ID_CODE) %>%
    mutate(cpralag = lag(CPRA)) %>%
    filter(cpralag != CPRA | is.na(cpralag)) %>% 
    select(-cpralag, -CHG_TIME) %>%
    ungroup()
  
  merged2 <- merged2 %>%
    inner_join(cpraWL, by = c("WL_ID_CODE")) %>%
    group_by(uniqueID) %>%
    mutate(innerrow = row_number()) %>%
    arrange(uniqueID, desc(innerrow)) %>%
    summarize(cpra = CPRA[predictDate >= CHG_DATE][1])
  
  ##################### 
  # Compute previously active
  #####################
  
  toMerge <- KidneyWaitlist %>% 
    select(WL_ID_CODE, CHG_DATE, nextDate, 
           UNOS_CAND_STAT_CD, totalDays)
  
  merged3 <- outputSet %>% select(-lastDate, -firstDate, -num3months) %>%
    mutate(uniqueID = row_number()) %>% tbl_df %>%
    inner_join(toMerge) %>%
    filter(predictDate >= CHG_DATE) %>%
    mutate(startWithin = predictDate <= nextDate)
  
  merged3$nextDate[merged3$startWithin] = merged3$predictDate[merged3$startWithin]
  
  merged3 <- merged3 %>% 
    mutate(totalDays = ifelse(startWithin,
                              as.numeric(difftime(nextDate, CHG_DATE, units = "days")),
                              totalDays)) %>%
    group_by(uniqueID) %>%
    summarize(totalActivePrev = sum(totalDays[UNOS_CAND_STAT_CD == 4010]))
  
  ##################### 
  # Compute future active for predRangeAmount
  #####################
  
  future <- merged4 %>% 
    inner_join(toMerge) %>% tbl_df() %>%
    filter(predictDate <= nextDate) %>% #remove everything that came before
    mutate(endDate = predictDate + days(predRangeAmount), 
           CHGBeyond = CHG_DATE > endDate,
           CHGFirst = CHG_DATE < predictDate,
           nextAfter = nextDate > endDate)
  
  future$CHG_DATE[future$CHGFirst] = future$predictDate[future$CHGFirst]
  
  future <- future %>%
    filter(!CHGBeyond)
  future$nextDate[future$nextAfter] = future$endDate[future$nextAfter]
  future <- future %>% 
    filter(CHG_DATE < nextDate) %>%
    mutate(totalDays = as.numeric(difftime(nextDate, CHG_DATE, units = "days")),
           totalDays) %>%
    group_by(uniqueID) %>%
    summarize(totalActiveNext = sum(totalDays[UNOS_CAND_STAT_CD == 4010]))
  
  
  ##################### 
  # Merge all data frames
  #####################
  
  combined <- merged %>% inner_join(original) %>%
    inner_join(merged2) %>%
    inner_join(merged3) %>%
    inner_join(future) %>%
    inner_join(allPatients %>% select(-numAccepted,
                                      -dr1, -dr2,
                                      -numOffers,
                                      -lastmatch))
  
  ##################### 
  # Begin preparing data frame for random forest models
  #####################
  
  TRAINING_START_DATE
  TRAINING_END_DATE
  VALIDATION_START_DATE
  VALIDATION_END_DATE
  TESTING_START_DATE
  TESTING_END_DATE
  predRangeAmount #in days
  OPOs
  ABOs
  
  combinedCopy <- combined
  
  #save.image(sprintf("preProcess-%d.Rdata", timeIndex))
  
  # prepare the data frame for predictive models
  combined <- combinedCopy
  combined <- combined %>%
    select(-uniqueID) %>%
    filter(
      totalActivePrev > predRangeAmount,
      (totalActiveNext > predRangeAmount * .9) | 
        (acceptsNext == 1 & as.numeric(daysUntilDesired) <= predRangeAmount),
      initdate <= predictDate,
      initdate >= as.POSIXct("2005-05-01", format = "%Y-%m-%d"),
      abo %in% c("O", "A", "AB", "B")) %>%
    mutate(receivedDesiredQ = factor(as.numeric(daysUntilDesired) <= predRangeAmount),
           years = as.numeric(difftime(predictDate, initdate, units = "days"))/365,
           q12intensity = numq12seen / totalActivePrev,
           medIntensity = numMediumSeen / totalActivePrev,
           lowIntensity = numLowSeen / totalActivePrev,
           curr_age = age + as.numeric(difftime(predictDate,
                                                initdate,
                                                units = "weeks"))/52,
           adult = (curr_age > 35 | is.na(curr_age)),
           pediatric = !adult,
           month = factor(month(predictDate)),
           cpraFactor = as.numeric(cpra < 0.8),
           avgOSQprev = ifelse(is.na(avgOSQprev),
                               500,
                               avgOSQprev),
           avgKDPIprev = ifelse(is.na(avgKDPIprev),
                                1,
                                avgKDPIprev),
           avgDRMISprev = ifelse(is.na(avgDRMISprev),
                                 2,
                                 avgDRMISprev)) %>%
    select(contains("intensity"),
           contains("avg"),
           cpra,
           cpraFactor,
           years,
           pediatric,
           receivedDesiredQ,
           predictDate,
           month,
           opo,
           abo)
  
  combined$opo <- factor(combined$opo)
  combined$abo <- factor(combined$abo)
  
  # divide into training, validation, testing sets
  tr <- combined %>% 
    filter(predictDate >= TRAINING_START_DATE,
           predictDate <= TRAINING_END_DATE) %>%
    select(-predictDate)
  va <- combined %>%
    filter(predictDate >= VALIDATION_START_DATE,
           predictDate <= VALIDATION_END_DATE) %>%
    select(-predictDate)
  ts <- combined %>%
    filter(predictDate >= TESTING_START_DATE,
           predictDate <= TESTING_END_DATE) %>%
    select(-predictDate)
  
  # remove the months column if we do not have
  # enough data to predict as factor
  levels.ts <- levels(factor(ts$month))
  levels.va <- levels(factor(va$month))
  levels.tr <- levels(factor(tr$month))
  
  
  if (length(setdiff(levels.ts, levels.tr)) != 0 |
      length(setdiff(levels.va, levels.tr)) != 0)
  {
    tr$month <- NULL
    va$month <- NULL
    ts$month <- NULL
  }
  
  # check boundary cases
  if(length(table(factor(tr$receivedDesiredQ))) == 1 |
     length(table(factor(va$receivedDesiredQ))) == 1 |
     length(table(factor(ts$receivedDesiredQ))) == 1 |
     nrow(tr) == 0 | 
     nrow(va) == 0 | 
     nrow(ts) == 0) {
    df = rbind(df, 
               c(timeIndex, predRangeAmount,
                 "all", "all",
                 nrow(tr) + nrow(va), sum(tr$receivedDesiredQ == TRUE) + sum(va$receivedDesiredQ == TRUE),
                 1, 
                 nrow(ts), sum(ts$receivedDesiredQ == TRUE),
                 NA, NA, NA, NA, NA))
    next
  }
  
  #######################################################
  
  aucList = list()
  nRows = list()
  
  
  nRows[[1]] = nrow(tr) + nrow(va)
  nRows[[2]] = sum(rbind(tr,va)$receivedDesiredQ == TRUE)
  nRows[[3]] = nrow(ts)
  nRows[[4]] = sum(ts$receivedDesiredQ == TRUE)
  
  ###################
  # logistic regression
  ###################
  log.mod <- glm(receivedDesiredQ ~., data = tr,
                 family = "binomial")
  pred.log <- predict(log.mod, newdata = va, type = "response")
  log.ROC = prediction(pred.log, va$receivedDesiredQ)
  valid.log.auc <- as.numeric(performance(log.ROC, "auc")@y.values)
  
  test.log.mod <- glm(receivedDesiredQ ~., data = rbind(tr, va),
                      family = "binomial")
  pred.log <- predict(log.mod, newdata = ts, type = "response")
  log.ROC = prediction(pred.log, ts$receivedDesiredQ)
  test.log.auc <- as.numeric(performance(log.ROC, "auc")@y.values)
  
  aucList["valid.log.mod"] = valid.log.auc
  aucList["test.log.mod"] = test.log.auc
  ###################
  
  ###################
  # classification trees
  ###################
  valid.auc = 0.5
  incumbent = 0
  for (cpValue in seq(0, 0.01, 0.001))
  {
    cart.mod <- rpart(receivedDesiredQ ~ ., data = tr, cp = cpValue)
    cart.preds <- predict(cart.mod, newdata = va, type = "prob")[, 2]
    ROCpred = prediction(cart.preds, va$receivedDesiredQ)
    check.auc <- as.numeric(performance(ROCpred, "auc")@y.values)
    if (check.auc > valid.auc) {
      incumbent = cpValue
      valid.auc = check.auc
    }
  }
  
  cart.mod <- rpart(receivedDesiredQ ~ ., data = rbind(tr, va), cp= incumbent)
  cart.preds <- predict(cart.mod, newdata = ts, type = "prob")[, 2]
  ROCpred = prediction(cart.preds, ts$receivedDesiredQ)
  cart.auc <- as.numeric(performance(ROCpred, "auc")@y.values)
  
  aucList["cart.valid.mod"] = valid.auc
  aucList["cart.test.mod"] = cart.auc
  ###################
  
  ###################
  # random forest
  ###################
  validateParam <- validateRF(tr, va)
  AggTR = rbind(tr, va)
  rf.valid <- randomForest(receivedDesiredQ ~ ., data = tr, 
                           mtry = validateParam[2],
                           ntree = validateParam[3])
  rf.preds <- predict(rf.valid, newdata = va, type = "prob")[, 2]
  
  if (all(rf.preds == rf.preds[1])){
    valid.auc <- 0
  } else {
    ROCpred = prediction(rf.preds, va$receivedDesiredQ)
    valid.auc <- as.numeric(performance(ROCpred, "auc")@y.values)
  }
  
  aucList["valid.rf.mod"] = valid.auc
  
  rf.mod <- randomForest(receivedDesiredQ ~ ., data = AggTR,
                         mtry = validateParam[2],
                         ntree = validateParam[3])
  rf.test <- predict(rf.mod, newdata = ts, type = "prob")[, 2]
  if (all(rf.test == rf.test[1])){
    test.auc <- 0
  } else {
    ROCpred = prediction(rf.test, ts$receivedDesiredQ)
    test.auc <- as.numeric(performance(ROCpred, "auc")@y.values)
  }
  
  aucList["test.rf.mod"] = test.auc
  
  save(rf.mod, file=sprintf("Models/agg-rf-%d.Rdata", timeIndex))
  
  varImpPlot(rf.mod)
  
  x <- rf.mod$importance
  rf.varImp.row <- t(x)[, 1:nrow(x)]
  rf.varImp.row <- data.frame(t(rf.varImp.row))
  
  opo = "all"
  abo = "all"
  predictOPO = "all"
  predictABO = "all"
  timingInd = timeIndex
  
  rf.varImp.row <- cbind(rf.varImp.row, data.frame(timingInd))
  
  df.var.Imp <- bind_rows(df.var.Imp, rf.varImp.row)
  
  ###################
  
  
  
  #######################################################
  
  df = rbind(df, 
             c(timeIndex, predRangeAmount,
               predictOPO, predictABO,
               nRows[[1]], nRows[[2]],
               nRows[[3]], nRows[[4]],
               aucList[[1]], aucList[[2]], aucList[[3]], aucList[[4]], aucList[[5]], aucList[[6]]))
}
  
save(df, file = sprintf("agg-df-%s", ID))
save(df.var.Imp, file = sprintf("agg-var-%s", ID))



