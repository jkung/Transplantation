##################### 
# Create the total amount of data that we need for
# each prediction time interval, plus the amount
# of data that we give to training, validation, and testing
#####################

#####################
# generate a data frame where each row
# is a configuration of date ranges to predict over
#####################
predictionRanges = c(91, 182, 365)
trainingData = c(365, 365 + 182, 
                 365 + 365, 365 + 365 + 182, 
                 365 + 365 + 365, 365 + 365 + 365 + 182)
trainingMonths = c(12, 18, 24, 30, 36, 42)
validationData = c(182, 365)
validationMonths = c(6, 12)
testingData = c(182, 365)
testingMonths = c(6, 12)

predRange = rep(predictionRanges, each = length(trainingData) * length(validationData) * length(testingData))
trainData = rep(trainingData,     each = length(validationData) * length(testingData), times = length(predictionRanges))
trainMonths = rep(trainingMonths, each = length(validationData) * length(testingData), times = length(predictionRanges))
validData = rep(validationData,     each = length(testingData), times = length(predictionRanges) * length(trainingData))
validMonths = rep(validationMonths, each = length(testingData), times = length(predictionRanges) * length(trainingData))
testData = rep(testingData,     times = length(predictionRanges) * length(trainingData) * length(validationData))
testMonths = rep(testingMonths, times = length(predictionRanges) * length(trainingData) * length(validationData))
#####################


#####################
# filter out illogical configurations (e.g. longer than the entire dataset,
# validation data exceeds testing data) and then rename columns
# to something more intuitive
#####################
timings = data.frame(predRange, trainData, validData, testData, trainMonths, validMonths, testMonths) %>%
  mutate(totalData = predRange * 3 + trainData + validData + testData,
         totalYears = totalData / 365) %>%
  filter(totalYears <= 6.1,
         testData >= validData) %>%
  mutate(predMonths = ifelse(predRange == 91, 3,
                             ifelse(predRange == 182, 6,
                                    ifelse(predRange == 365, 12, NA)))) %>%
  mutate(monthsToYears = (predMonths*3 + trainMonths + validMonths + testMonths)/ 12) %>%
  arrange(totalData)
#####################

# the dates tested in the first part of the paper correspond to indices 32, 34, 38