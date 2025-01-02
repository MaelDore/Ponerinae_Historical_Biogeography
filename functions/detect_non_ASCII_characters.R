####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

## Goal

# Detect all unique non-ASCII characters in a data.frame and print them

## Input

# data.frame

## Output

# Print the list of unique non-ASCII characters found in the data.frame


## Detect non-ASCII characters in a data.frame

detect_non_ASCII_characters_from_df <- function (df)
{
  all_chars <- paste(unlist(df), collapse = "")
  all_non_ASCII_chars <- unique(unlist(str_extract_all(string = all_chars, pattern = "[^\t\n\r\x20-\x7E]")))
  cat(paste0(length(all_non_ASCII_chars), " unique non-ASCII characters detected: ", paste(all_non_ASCII_chars, collapse = " ")))
}