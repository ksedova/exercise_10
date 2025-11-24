rm(list=ls())

setwd('D:/university/Amasters/firstsemester/PRG/exercise_10')

# Suffix Trees

### Task 1
# * In R, implement a function `SuffixArray()` to create a suffix array from a string.
# 
# * Input:
#   * A `DNAString` object .
# 
# * Output:
#   * A vector of integers.

library(Biostrings)

SuffixArray <- function(dna_string) {
  # 1. Validate input is a DNAString
  if (!is(dna_string, "DNAString")) {
    stop("Input must be a DNAString object.")
  }
  
  # 2. Convert to standard character for manipulation
  # (Note: For massive genomes, specific C-level libraries are preferred 
  # to avoid memory overflow, but this is the standard R implementation logic)
  s_char <- as.character(dna_string)
  n <- nchar(s_char)
  
  # 3. Generate all suffixes
  # substring(x, i) extracts from index i to the end of the string
  suffixes <- substring(s_char, 1:n)
  
  # 4. Return the order (this creates the Suffix Array)
  # 'order' returns the indices that would sort the vector
  return(order(suffixes))
}

# Create a DNAString object
dna_input <- DNAString("ACGAAC")

# Run the function
sa_result <- SuffixArray(dna_input)

# Print the result
print(sa_result)


### Task 2
# * In R, implement a function `InverseSuffixArray()` to create an inverse suffix array from a suffix array.
# 
# * Input:
#   * A vector of integers representing suffix array.
# 
# * Output:
#   * A vector of integers.

InverseSuffixArray <- function(suffix_array) {
  # 1. Validate input
  if (!is.integer(suffix_array) && !is.numeric(suffix_array)) {
    stop("Input must be a vector of integers.")
  }
  
  n <- length(suffix_array)
  
  # 2. Pre-allocate the vector for efficiency
  isa <- integer(n)
  print(isa)
  
  # 3. Compute Inverse Suffix Array
  # Since R supports vectorized indexing, we can assign ranks in one step.
  # We map the index (1:n) to the positions defined by the suffix_array.
  isa[suffix_array] <- 1:n
  
  return(isa)
}

# Run the function
isa_result <- InverseSuffixArray(sa_result)

# Print the result
print(isa_result)

### Task 3
# * In R, implement a function `LCPArray()` according to pseudocode.
# 
# * Input:
#   * `text` A `DNAString` representing analyzed string.
# * `SA` A vector of integers representing a suffix array.
# * `ISA` A vector of integers representing an inverse suffix array.
# 
# * Output:
#   * `LCP` A vector of integers.
# 
# ```
# LCPArray(text, SA, ISA)
# 1   LCP[1] <- -1
# 2   LCP[m + 1] <- -1
# 3   l <- 0
# 4   for i <- 1 to m
# 5     j <- ISA[i]
# 6     if j > 1 then
# 7       k <- SA[j - 1]
# 8       while text[k + l] = text[i + l] do
# 9         l <- l + 1
# 10      LCP[j] <- l
# 11      l <- max{l - 1, 0}
# 12  return LCP
# ```
# 
# **Hint:** 
#   The text will be indexed at *m* + 1 position, that does not exist. Add one character at the end of the text
# (in general use `$`, for `DNAString` in R use `+`).
# 

LCPArray <- function(text, SA, ISA){
  m <- length(text)
  
  # --- Handling the Hint ---
  # Convert DNAString to a standard character vector for easy indexing 
  # and append the sentinel "+" to handle the boundary check in the while loop.
  s_vec <- strsplit(as.character(text), "")[[1]]
  s_vec <- c(s_vec, "+")
  
  LCP <- integer(m + 1)
  
  LCP[1] <- -1
  LCP[m + 1] <- -1
  
  l <- 0
  
  for (i in 1:m){
    j <- ISA[i]
    if (j > 1){
      k <- SA[j - 1]
      while (s_vec[k + l] == s_vec[i + l]){
        l <- l + 1
      }
      LCP[j] <- l
      l <- max(l - 1, 0)
    }
  }
  return(LCP)
}

lcp <- LCPArray(dna_input, sa_result, isa_result)
print(paste("LCP:", paste(lcp, collapse=" ")))

### Task 4
# * In R, implement a function `BinarySearchSA()` according to the following pseudocode.
# 
# * Input:
#   * `pattern` A `DNAString` representing a pattern to be found.
# * `text` A `DNAString` representing a text to be searched.
# * `SA` A vector of integers representing a suffix array of `text`.
# 
# * Output:
#   * A vector of two integers (the first and the last indexes of suffix array, where the pattern was found).
# 
# ```
# BinarySearchSA(pattern, text, SA)
# 1   minIndex <- 1
# 2   maxIndex <- length (text)
# 3   while minIndex < maxIndex
# 4     midlIndex <- floor(minIndex + maxIndex) / 2
# 5     if pattern <= suffix of text starting at position SA(midlIndex)
# 6       maxIndex <- midlIndex
# 7     else
#   8       minIndex <- midlIndex + 1
# 9   First <- minIndex
# 10  maxIndex <- length(text)
# 11  while maxIndex > minIndex
# 12    midlIndex <- floor(minIndex + maxIndex) / 2
# 13    if suffix of text starting at position SA(midlIndex) <= pattern
# 14      minIndex <- midlIndex + 1
# 15    else
#   16      maxInd <- midlIndex
# 17  Last <- maxIndex - 1
# 18  if Last < First
# 19    return('Pattern does not appear in text')
# 20  else
#   21    return First, Last
# ```

BinarySearchSA <- function(pattern, text, SA){
  # 1. Setup
  # Ensure inputs are proper types
  if (is(pattern, "DNAString")) pattern <- as.character(pattern)
  if (is(text, "DNAString")) text <- as.character(text) else text <- text
  
  p_len <- nchar(pattern)
  
  minIndex <- 1
  maxIndex <- length(SA)
  
  while(minIndex < maxIndex){
    midlIndex <- floor((minIndex + maxIndex) / 2)
    
    # 1. Look at the Suffix Array to find where the suffix starts in the text
    suffix_start_pos <- SA[midlIndex]
    
    # 2. Extract the actual DNA letters from the text
    # We only take the first P_len characters to compare prefixes
    suffix_string <- substring(text, suffix_start_pos, suffix_start_pos + p_len - 1)
    
    if (pattern <= suffix_string){
      maxIndex <- midlIndex
    }else{
      minIndex <- midlIndex + 1
    }
  }
  First <- minIndex
  maxIndex <- length(SA)
  
  while (maxIndex > minIndex){
    midlIndex <- floor((minIndex + maxIndex) / 2)
    
    # Extract suffix again for comparison
    suffix_start_pos <- SA[midlIndex]
    suffix_string <- substring(text, suffix_start_pos, suffix_start_pos + p_len - 1)
    
    if (suffix_string <= pattern){
      minIndex <- midlIndex + 1
    }else{
      maxIndex <- midlIndex
    }
  }
  Last <- maxIndex - 1
  if (Last < First){
    return('Pattern does not appear in text')
  }else{
    return(c(First, Last))
  }
}

# 1. Setup Data
dna_input <- DNAString("CTAATAATG")
# Recall SA for "ACGAAC" is: 4 5 1 6 2 3
sa <- SuffixArray(dna_input)

# 2. Define a pattern to search
# "AC" appears twice: at index 1 ("ACGAAC") and index 5 ("AC")
pattern <- DNAString("GA")

# 3. Run Search
result <- BinarySearchSA(pattern, dna_input, sa)

# 4. Print Result
print(result)