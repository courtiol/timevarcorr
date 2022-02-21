## code to prepare `stockprice` dataset goes here

## data downloaded from: https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices-0.2-final.zip

stockprice_raw <- readr::read_delim(file = "./data-raw/OxfordManRealizedVolatilityIndices.csv", skip = 3)
stockprice <- stockprice_raw[, c("DateID", "SPX2.r", "FTSE2.r", "N2252.r", "GDAXI2.r", "IXIC2.r")]
stockprice$DateID <- lubridate::ymd(stockprice$DateID)
colnames(stockprice)[colnames(stockprice) == "SPX2.r"] <- "SP500"
colnames(stockprice)[colnames(stockprice) == "FTSE2.r"] <- "FTSE100"
colnames(stockprice)[colnames(stockprice) == "N2252.r"] <- "Nikkei"
colnames(stockprice)[colnames(stockprice) == "GDAXI2.r"] <- "DAX"
colnames(stockprice)[colnames(stockprice) == "IXIC2.r"] <- "NASDAQ"
stockprice <- stockprice[stockprice$DateID >= "2000-04-01" & stockprice$DateID <= "2017-12-30", ]
stockprice$Event <- with(stockprice, dplyr::case_when(DateID >= "2001-04-01" & DateID <= "2001-12-27" ~ "dot com bubble",
                                                      DateID >= "2007-11-21" & DateID <= "2009-05-22" ~ "global financial crisis",
                                                      DateID >= "2011-06-05" & DateID <= "2012-05-30" ~ "European debt crisis",
                                                      TRUE ~ "baseline"))

usethis::use_data(stockprice, overwrite = TRUE)

#plot(stockprice$SP500 ~ stockprice$NASDAQ)
