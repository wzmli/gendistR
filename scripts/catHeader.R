# Function to create dynamic tabsets:

catHeader <- function(text = "", level = 3) {
  cat(paste0("\n\n", 
             paste(rep("#", level), collapse = ""), 
             " ", text, "\n"))
}

#See: https://stackoverflow.com/questions/53444333/dynamic-tabsets-with-multiple-plots-r-markdown/53444654#53444654

# Nathan K Smith
# March 2023