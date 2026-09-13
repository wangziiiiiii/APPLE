# Names evaluated within dplyr/data.table data masks, not package globals.
utils::globalVariables(c(
  ".", "cluster_id", "condition", "five_prime_end", "lib_id", "p_value",
  "pac_id", "q_value", "tail_length", "tail_lengths", "test_length"
))
