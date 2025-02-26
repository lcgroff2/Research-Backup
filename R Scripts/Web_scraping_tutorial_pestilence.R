#Install Selenium from CRAN
install.packages("RSelenium")
#Install web driver manager from CRAN 
install.packages("wdman")

#Libraries
library(RSelenium)
library(wdman)

#Open a chrome browser named remDr
driver<- rsDriver(browser=c("chrome"),chromever = "80.0.3987.16", port = 4471L)
remDr <- driver[["client"]]
remDr$maxWindowSize()

#Let's go to the NC dept of health
url <- "https://www.ncdhhs.gov/divisions/public-health/covid19/covid-19-nc-case-count"
remDr$navigate(url)

#We want to find what information is available by county and click on it
webElem <- remDr$findElement(using = "id", value = 'by-counties')
webElem$clickElement()

#Here we find and grab the table showing the progression of the pestilence
webElem <- remDr$findElement(using= "xpath", 
                             value='//*[@id="ui-accordion-ui-id-1-panel-2"]/section')
elemtxt <- webElem$getElementText()

#Split the data using end of line \n
split_text <- strsplit(elemtxt[[1]],"\n")
#Remove data we accidentally grabbed earlier
split_text2 <- split_text[[1]][3:(length(split_text[[1]])-1)]
#Initialize a dataframe 
Happy_dataframe <- data.frame(County = 1:length(split_text2))
#Loop over all counties and place data in appropriate places
for (Jon in 1:length(split_text2)) {
  Happy_dataframe[Jon,"County"] <- 
    strsplit(split_text2[Jon]," ")[[1]][1]
  Happy_dataframe[Jon,"Laboratory-Confirmed Cases"] <- 
    strsplit(split_text2[Jon]," ")[[1]][3]
  Happy_dataframe[Jon,"Deaths"] <- 
    strsplit(split_text2[Jon]," ")[[1]][4]
}

View(Happy_dataframe)










