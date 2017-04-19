# re-order DR antigens from small to large
DRmatchingSorted<- DRmatching %>% mutate(sDDR = pmin(DDR1, DDR2),
                                    lDDR = pmax(DDR1, DDR2),
                                    sDR = pmin(DR1, DR2),
                                    lDR = pmax(DR1, DR2),
                                    DDR1 = sDDR,
                                    DDR2 = lDDR,
                                    DR1 = sDR,
                                    DR2 = lDR,
                                    DR1 = ifelse(DR1 == 0, NA, DR1),
                                    DR2 = ifelse(DR2 == 0, NA, DR2),
                                    DDR1 = ifelse(DDR1 == 0, NA, DDR1),
                                    DDR2 = ifelse(DDR2 == 0, NA, DDR2)) %>% 
  select(-sDDR, -lDDR, -sDR, -lDR)

# create data frame with just donors that 
# were observed in the training set
matching <- donors %>% filter(date < TRAINING_END_DATE) %>% 
  select(DONOR_ID, opo) %>% 
  inner_join(DRmatchingSorted %>% select(DONOR_ID, DDR1, DDR2)) %>%
  group_by(DONOR_ID) %>% summarize(opo = first(opo),
                                   DDR1 = first(DDR1),
                                   DDR2 = first(DDR2)) %>%
  mutate(DDR2 = ifelse(DDR2 == 103, 1, DDR2),
         DDR1 = ifelse(DDR1 == 103, 1, DDR1)) %>% 
  group_by(opo, DDR1, DDR2) %>%
  summarize(counts = n()) %>%
  filter(!is.na(DDR1),
         !is.na(DDR2))

# compute average mismatch for each (DR1, DR2, opo) combination
temp <- allPatientsPre %>% group_by(dr1, dr2, opo) %>% summarize(npatients = n()) %>%
  filter(!is.na(dr1), !is.na(dr2), !is.na(opo))
out <- mapply(getAvgMismatch, temp$dr1, temp$dr2, temp$opo)
temp$avgMIS <- out
rm(out)

allPatients <- allPatientsPre %>% inner_join(temp %>% select(dr1, dr2, opo, avgMIS))






