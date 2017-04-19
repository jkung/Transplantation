#####################
# Load all relevant libraries
#####################
library(dplyr)
library(lubridate)
library(rpart)
library(rpart.plot)
library(ROCR)
library(randomForest)
library(beepr)
library(ggplot2)
library(caret)
library(e1071)
library(extrafont)

#####################
# function getPatientsKidneys():
# returns the patients and kidneys data frames in a list
# with pertinent information summarized
#####################
getPatientsKidneys<- function(){
  patients <- PTRdata %>% group_by(WL_ID_CODE) %>% 
    summarize(initdate = first(INIT_DATE),
              numAccepted = sum(ACCEPTED),
              opo = first(LISTING_CTR_DSA),
              age = first(AGE),
              dr1 = first(DR1),
              dr2 = first(DR2),
              numOffers = n(),
              abo = first(ABO_CAND),
              lastmatch = max(MATCH_DATE)) %>%
    arrange(initdate) 
  
  kidneys <- PTRdata %>% group_by(DONOR_ID) %>% 
    summarize(date = first(MATCH_DATE),
              numAccepted = sum(ACCEPTED),
              opo = first(OPO_CTR),
              age = first(AGE_DON),
              kdpi = first(KDPI_CALC),
              maxOffers = max(ORGAN_SEQUENCE_NUM),
              ecd = first(ECD_DONOR),
              abo = first(ABO_DON),
              quality = first(kidneyquality)) %>%
    arrange(date)
  return(list(patients, kidneys))
}

#####################
# function addNextDate():
# takes removeRedundantRows and adds the nextDate column
# also removes all entries before 2007-05-01 or if an
# entry spans that date, makes the start date 2007-05-01
#####################
addNextDate <- function(){
  copy <- KidneyWaitlist
  
  copy$CHG_DATE <- mdy(copy$CHG_DATE)
  copy <- copy %>% group_by(WL_ID_CODE) %>% 
    mutate(sameCD = (UNOS_CAND_STAT_CD == lag(UNOS_CAND_STAT_CD))) %>%
    filter(!sameCD | is.na(sameCD) | CHG_TY == "D") %>% 
    mutate(nextDate = lead(CHG_DATE))

  a <- as.POSIXct("2013-06-30", format = "%Y-%m-%d")
  
  test <- copy$nextDate
  test2 <- copy$CHG_TY
  test3 <- copy$CHG_DATE
  
  test[is.na(test) & test2 == "D"] = test3[is.na(test) & test2 == "D"]
  test[is.na(test) & (test2 == "M" | test2 == "A")] = a
  
  copy$nextDate = test
  
  return(copy)
}

#####################
# function filterDates(beginDate, endDate):
# takes a date range, and cuts off the statuses to be within these date ranges
# this is useful, for example, for different training and testing sets
#####################
filterDates <- function(beginDate, endDate){
  beginDate = as.POSIXct(beginDate, format = "%Y-%m-%d")
  endDate = as.POSIXct(endDate, format = "%Y-%m-%d")
  
  KidneyWaitlist <- KidneyWaitlist %>% filter(nextDate >= beginDate) 
  test <- KidneyWaitlist$CHG_DATE 
  test[test < beginDate] = beginDate
  KidneyWaitlist$CHG_DATE = test
  
  KidneyWaitlist <- KidneyWaitlist %>% filter(CHG_DATE <= endDate)
  test <- KidneyWaitlist$nextDate
  test[test > endDate] = endDate
  KidneyWaitlist$nextDate = test
  
  return(KidneyWaitlist)
}

#####################
# function validateRF(training set, validation set)
# builds a random forest model on the training set for
# parameters in the set:
#  mtry \in {3, .., nvar - 1}
#  ntree \in {500, 1000}
# and returns configurations that perform best on validation
#####################
validateRF <- function(tr, va){
  incumbent = c(0, NA, NA)
  total = (tr %>% names %>% length)-1
  for (mtry in seq(3,total))
    for (ntree in c(500, 1000))
    {
      rf.mod = randomForest(receivedDesiredQ ~., data = tr, mtry = mtry, ntree = ntree)
      rf.preds <- predict(rf.mod, newdata = va, type = "prob")[, 2]
      if (all(rf.preds == rf.preds[1]))
        next
      ROCpred = prediction(rf.preds, va$receivedDesiredQ)
      rf.auc <- as.numeric(performance(ROCpred,"auc")@y.values)
      if (rf.auc > incumbent[1]){
        incumbent = c(rf.auc, mtry, ntree)
      }  
    }
  
  return(incumbent)
}


#####################
# function numMismatch(recip DR1, recip DR2, donor DR1, donor DR2)
# that takes patient DR antigens, donor DR antigens (respectively)
# and returns the number of mismatches
#####################
numMismatch <- function(DR1, DR2, DDR1, DDR2){
  if (DR1 == 103) DR1 = 1
  if (DR2 == 103) DR2 = 1
  if(DDR1 == 103) DDR1 = 1
  if(DDR2 == 103) DDR2 = 1
  
  pDR <- c(DR1, DR2)
  dDR <- c(DDR1, DDR2)
  if (setdiff(pDR, dDR) %>% length == 0)
  {
    return(0)
  } else if (DDR1 %in% c(2, 3, 16, 17, 18, 6, 14, 1403, 1404) | 
             DDR2 %in% c(2, 3, 16, 17, 18, 6, 14, 1403, 1404))
  {
    added = 0
    if (DR1 == 14 | DR2 == 14){
      before = length(unique(pDR))
      pDR <- c(pDR, 14, 1403, 1404)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    if(DR1 == 16 | DR2 == 16){
      before = length(unique(pDR))
      pDR <- c(pDR, 16, 2)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    if(DR1 == 17 | DR2 == 17){
      before = length(unique(pDR))
      pDR <- c(pDR, 17, 3)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    if(DR1 == 18 | DR2 == 18){
      before = length(unique(pDR))
      pDR <- c(pDR, 18, 3)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    if(DR1 == 1403 | DR2 == 1403){
      before = length(unique(pDR))
      pDR <- c(pDR, 1403, 14, 6)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    if(DR1 == 1404 | DR2 == 1404){
      before = length(unique(pDR))
      pDR <- c(pDR, 1404, 14, 6)
      after = length(unique(pDR))
      added = added + after - before
    }
    
    return(length(setdiff(pDR, dDR)) - added)
  }
  else 
  {
    return(setdiff(dDR, pDR) %>% length)
  }
}

#####################
# function getAvgMismatch(recip DR1, recip DR2, opo)
# returns the average mismatch observed *in the training set*
# note that the data frame matching has already been filtered 
# to reflect only those observations in the training set
#####################
getAvgMismatch <- function(DR1, DR2, opoinput){
  opoinput <- as.character(opoinput)
  temp <- matching %>% ungroup() %>% mutate(opo = as.character(opo)) %>% filter(opo == opoinput)
  temp$DR1 <- DR1
  temp$DR2 <- DR2
  temp$DRMIS <- mapply(numMismatch, temp$DR1, temp$DR2, temp$DDR1, temp$DDR2)
  temp <- temp %>% mutate(numerator = counts * DRMIS)
  numerator <- sum(temp$numerator)
  denom <- sum(temp$counts)
  return(numerator/denom)
}

