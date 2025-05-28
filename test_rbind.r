my_df <- data.frame(A = 1:3, B = letters[1:3])
# Trace memory of my_df
tracemem(my_df)
# Append a new row using rbind()
my_df <- rbind(my_df, data.frame(A = 4, B = "d"))

# tracemem[<old_address> -> <new_address>] indicates that
# the original memory address of my_df was copied
# to a new location when rbind() created the larger data frame.
print(my_df)

# Append another row
my_df <- rbind(my_df, data.frame(A = 5, B = "e"))


library(tibble)
library(dplyr)
my_df <- tibble(A = 1:3, B = letters[1:3])
tracemem(my_df)
my_df <- bind_rows(my_df, tibble(A = 4, B = "d"))
tracemem(my_df)
print(my_df)
