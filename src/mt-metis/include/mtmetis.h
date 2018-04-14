/**
 * @file mtmetis.h
 * @brief Library entry points
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-01
 */


#ifndef MTMETIS_H
#define MTMETIS_H

#define MTMETIS_VER_MAJOR 0
#define MTMETIS_VER_MINOR 2
#define MTMETIS_VER_SUBMINOR 0

int mtmetis_partkway(const idx_t * nvtxs, const idx_t * xadj, 
    const idx_t * adjncy, const idx_t * vwgt, const idx_t * adjwgt,
    const idx_t * nparts, idx_t * where, idx_t * edgecut);

#endif
