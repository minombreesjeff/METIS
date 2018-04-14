# {#mainpage}

mt-Metis
=============================

mt-Metis is a multithreaded multilevel graph partitioning an ordering tool. It
is based on the algorithms used in Metis and ParMetis 
(http://cs.umn.edu/~metis).


The 'mtmetis' Executable
-----------------------------

To partition a graph with mt-Metis into 16 parts:

    mtmetis test.graph 16 test.part

To generate a nested dissection ordering a of a graph/matrix:

    mtmetis -p nd test.graph test.perm  

Many other runtime parameters are available. Use 
the '-h' option to view these.

    mtmetis -h

Be warned that many of the options are for experimentation purposes and may
signficantly effect the performance and functionality of mt-Metis.
    


The mt-Metis API
-------------------------------

The file [mtmetis.h](@ref mtmetis.h) is the header that should be included
by external programs wishing link to mt-Metis. There are two high level
functions, mtmetis_partkway() for partitioning, and mtmetis_nd() for generating
orderings. At this time mt-Metis is highly experimental, and its
API is subject to change.

The following is an example of how to partition an unweighted graph stored in 
CSR format in the arrays xadj and
adjncy, with n vertices, and we want to generate k parts.

    mtmetis_partkway(n,xadj,adjncy,NULL,NULL,k,part,edgecut);

In this example, part is an array of length n that will have a partition ID
assigned to each vertex (i.e., vertex v is assigned to partition part[v]).


The following is an example of how to order an unweighted graph/matrix in 
similar CSR format (with n rows/columns).

    mtmetis_nd(n,xadj,adjncy,NULL,NULL,perm);

In this example, perm is an array of length n that will have the new row ID
and column ID assigned to each row/column (i.e., row r's index in the new
ordering is perm[r]).



These functions are documented [here](@ref mtmetis.h).

