## RSQLite safe join
.fullJoin <- function(table1, table2) {
	join1 <- left_join(table1, table2)
	join2 <- left_join(table2, table1)
	union(join1, join2)
}
