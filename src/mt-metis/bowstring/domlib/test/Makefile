CC = gcc
RM = rm -rfv

CSTD = -std=gnu99
ifdef DEBUG
    OPTOPTS = -g -DDEBUG
else
    OPTOPTS = -march=native -O3
endif

COPTS = $(OPTOPTS) $(CSTD) -Werror -Wmissing-prototypes \
        -Wall -Wextra -Wno-unused-parameter -Wno-unused-function \
        -Wmissing-declarations

LIBS = -L../ -ldom -fopenmp -lm
HDRS = -I../

BUFFER_TESTS = buffer_add.pass \
               buffer_clearreset.pass


CMDLINE_TESTS = cmdline_parse.pass

DLDPQ_TESTS = dldpq_fill.pass \
              dldpq_inc.pass \
              dldpq_remove.pass

FILE_TESTS = file_get_next_line.pass

HEAP_TESTS = heap_push.pass \
             heap_pop.pass \

PRIORITYQUEUE_TESTS = \
             priority_queue_pushpop.pass \
             priority_queue_update.pass \
             priority_queue_delete.pass \
             priority_queue_clear.pass \
             priority_queue_peek.pass \
             priority_queue_max.pass

LIST_TESTS = list_linked_create.pass \
             list_linked_pushpop.pass \
             list_linked_get.pass \
             list_linked_replace.pass

MACROS_TESTS = macros_clz.pass

MAP_TESTS = hashtable_create.pass \
            hashtable_get.pass \
            hashtable_put.pass \
            hashtable_remove.pass \
            hashtable_contains.pass

MATH_TESTS = math_incset.pass \
             math_cyclicperm.pass \
             math_blockcyclicperm.pass \
             math_prefixsum_inc.pass \
             math_prefixsum_exc.pass \
             math_prefixshift.pass \
             math_differentiate.pass \
             math_max_index.pass \
             math_max_value.pass \
             math_set_max.pass \
             math_min_index.pass \
             math_min_value.pass \
             math_set_min.pass \
             math_sum.pass \
             math_product.pass \
             math_downpow2.pass \
             math_uppow2.pass \
             math_stable_sum.pass \
             math_intersection_size.pass \
             math_reversebits.pass

MEM_TESTS = mem_calloc.pass \
            mem_set.pass \
            mem_copy.pass \
            mem_realloc.pass \
            mem_duplicate.pass \
            mem_r_sym_calloc.pass \
            mem_r_dalloc.pass

OMP_TESTS = omp_combuffer_create.pass \
            omp_combuffer_add.pass \
            omp_combuffer_send.pass \
            omp_combuffer_clear.pass \
            omp_reduction_init.pass \
            omp_reduction_sumareduce.pass

RAND_TESTS = rand_pseudo_shuffle.pass \
             rand_shuffle.pass

SEARCH_TESTS = search_binary.pass

SET_TESTS = set_iset_add.pass \
            set_iset_remove.pass \
            set_iset_remove_index.pass \
            set_iset_populate.pass \
            set_djset_findjoin.pass

SORT_TESTS = sort_counting_v.pass \
             sort_insertion.pass \
             sort_quick.pass \
             sort_radix.pass \
             sort_counting_v.pass 

STATS_TESTS = stats_arithmetic_mean.pass

STRING_TESTS = string_endswith.pass

TREE_TESTS = tree_addget.pass

UTIL_TESTS = util_bytes.pass

TESTS = $(MEM_TESTS) \
        $(MACROS_TESTS) \
        $(MATH_TESTS) \
        $(OMP_TESTS) \
        $(SEARCH_TESTS) \
        $(SORT_TESTS) \
        $(STATS_TESTS) \
        $(UTIL_TESTS) \
        $(SET_TESTS) \
        $(HEAP_TESTS) \
        $(PRIORITYQUEUE_TESTS) \
        $(CMDLINE_TESTS) \
        $(LIST_TESTS) \
        $(FILE_TESTS) \
        $(MAP_TESTS) \
        $(DLDPQ_TESTS) \
        $(RAND_TESTS) \
        $(STRING_TESTS) \
        $(TREE_TESTS)

all: $(TESTS)

test: $(TESTS)

# STUFF AND JUNK ##############################################################

%.pass: %.test
	@printf "Running $<... " && \
    ./$< && \
    touch $@ 

%.test: %.c test.h test.o ../*.h ../libdom.a
	$(CC) $(COPTS) $< test.o $(HDRS) $(LIBS) -o $@

test.o: test.h test.c ../*.h ../libdom.a
	$(CC) $(COPTS) -c test.c $(HDRS) $(LIBS) -o $@

clean:
	$(RM) *.pass
	$(RM) *.test
	$(RM) *.o
	$(RM) *.swp 

distclean: clean
	
rebuild: distclean all 

.PHONY: test all
