/**
 * @file GraphInHandle.hpp
 * @brief Class for reading all graph types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef WILDRIVER_GRAPHINHANDLE_HPP
#define WILDRIVER_GRAPHINHANDLE_HPP




#include <vector>
#include <memory>

#include "IGraphReader.hpp"




namespace WildRiver 
{


class GraphInHandle
{
  private:
    /**
     * @brief A pointer to the underlying graph reader.
     */
    std::shared_ptr<IGraphReader> reader;


    /**
     * @brief Private copy constructor declared to disable copying.
     *
     * @param handle The handle to copy.
     */
    GraphInHandle(
        GraphInHandle const & handle);


    /**
     * @brief Private assignment operator declared to disable copying.
     *
     * @param handle The handle to copy.
     *
     * @return The new handle.
     */
    GraphInHandle & operator=(
        GraphInHandle const & handle);


  public:
    /**
     * @brief Create a new file handle for reading matrices.
     *
     * @param fname The filename/path of the file to read.
     */
    GraphInHandle(
        std::string const & fname);


    /**
     * @brief Close the handle.
     */
    ~GraphInHandle();


    /**
     * @brief Read the CSR structure of the graph.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights (length nvtxs*nvwgt). This may be NULL in
     * order to ignore vertex weights. If it is specified and the file does not
     * contain vertex weights, it will be filled with ones.
     * @param adjwgt The edge weights (length nedges). This may be NULL in
     * order to ignore edge weights. If it is specified and the file does not
     * contain edge weights, it will be filled with ones.
     */
    void readGraph(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt);

  
    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts);


    /**
     * @brief Get the information of the next vertex.
     *
     * @param vwgt The vertex weight(s).
     * @param list The adjacency list of the vertex.
     *
     * @return The current vertex number.
     */
    dim_t getNextVertex(
        std::vector<val_t> & vwgt,
        std::vector<MatrixEntry> & list);


};




}




#endif
