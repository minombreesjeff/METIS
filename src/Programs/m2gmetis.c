/*
 * Copyright 1994-2011, Regents of the University of Minnesota
 *
 * m2gmetis.c
 *
 * Drivers for the mesh-to-graph coversion routines
 *
 * Started 5/28/11
 * George
 *
 * $Id: m2gmetis.c 10046 2011-06-01 14:13:40Z karypis $
 *
 */

#include "metisbin.h"



/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[])
{
  mesh_t *mesh;
  graph_t *graph;
  params_t *params;

  params = parse_cmdline(argc, argv);

  gk_startcputimer(params->iotimer);
  mesh = ReadMesh(params);

  gk_stopcputimer(params->iotimer);

  if (mesh->ncon > 1) {
    printf("*** Meshes with more than one balancing constraint are not supported yet.\n");
    exit(0);
  }

  M2GPrintInfo(params, mesh);

  graph = CreateGraph();

  gk_startcputimer(params->parttimer);
  switch (params->gtype) {
    case METIS_GTYPE_DUAL:
      METIS_MeshToDual(&mesh->ne, &mesh->nn, mesh->eptr, mesh->eind, 
            &params->ncommon, &params->numflag, &graph->xadj, &graph->adjncy);

      graph->nvtxs  = mesh->ne;
      graph->nedges = graph->xadj[graph->nvtxs]; 
      graph->ncon   = 1;
      break;

    case METIS_GTYPE_NODAL:
      METIS_MeshToNodal(&mesh->ne, &mesh->nn, mesh->eptr, mesh->eind, 
            &params->numflag, &graph->xadj, &graph->adjncy);

      graph->nvtxs  = mesh->nn;
      graph->nedges = graph->xadj[graph->nvtxs]; 
      graph->ncon   = 1;
      break;
  }
  gk_stopcputimer(params->parttimer);

  /* Write the graph */
  gk_startcputimer(params->iotimer);
  WriteGraph(graph, params->outfile);
  gk_stopcputimer(params->iotimer);

  M2GReportResults(params, mesh, graph);

  FreeGraph(&graph);
  FreeMesh(&mesh);
  gk_free((void **)&params->filename, &params->outfile, LTERM);
}


/*************************************************************************/
/*! This function prints run parameters */
/*************************************************************************/
void M2GPrintInfo(params_t *params, mesh_t *mesh)
{ 
  printf("******************************************************************************\n");
  printf("%s", METISTITLE);
  printf(" (HEAD: %s, Built on: %s, %s)\n", SVNINFO, __DATE__, __TIME__);
  printf(" size of idx_t: %zubits, real_t: %zubits, idx_t *: %zubits\n", 
      8*sizeof(idx_t), 8*sizeof(real_t), 8*sizeof(idx_t *));
  printf("\n");
  printf("Mesh Information ------------------------------------------------------------\n");
  printf(" Name: %s, #Elements: %"PRIDX", #Nodes: %"PRIDX"\n", 
      params->filename, mesh->ne, mesh->nn); 
  
  printf("Options ---------------------------------------------------------------------\n");
  printf(" gtype=%s, ncommon=%"PRIDX", outfile=%s\n", 
      gtypenames[params->gtype], params->ncommon, params->outfile);

  printf("\n");
}


/*************************************************************************/
/*! This function does any post-metis reporting */
/*************************************************************************/
void M2GReportResults(params_t *params, mesh_t *mesh, graph_t *graph)
{ 

  gk_startcputimer(params->reporttimer);

  printf(" - #nvtxs: %"PRIDX", #edges: %"PRIDX"\n", graph->nvtxs, graph->nedges);

  gk_stopcputimer(params->reporttimer);


  printf("\nTiming Information ----------------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3"PRREAL"\n", gk_getcputimer(params->iotimer));
  printf("  Partitioning: \t\t %7.3"PRREAL"   (METIS time)\n", gk_getcputimer(params->parttimer));
  printf("  Reporting:    \t\t %7.3"PRREAL"\n", gk_getcputimer(params->reporttimer));
  printf("******************************************************************************\n");

}
