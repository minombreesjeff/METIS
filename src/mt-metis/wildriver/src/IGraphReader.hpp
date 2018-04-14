/**
 * @file IGraphReader.hpp
 * @brief Interface for reading graphs.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef WILDRIVER_IGRAPHREADER_HPP
#define WILDRIVER_IGRAPHREADER_HPP





#include <string>
#include <vector>


#include "base.h"
#include "MatrixEntry.hpp"
#include "Exception.hpp"




namespace WildRiver
{


class IGraphReader
{
  private:


  protected:


  public:
    /**
     * @brief Destructor.
     */
    virtual ~IGraphReader() 
    {
    }


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
    virtual void read(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt) = 0;

  
    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    virtual void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts) = 0;


    /**
     * @brief Set the cursor to the first vertex in the file.
     */
    virtual void firstVertex() = 0;


    /**
     * @brief Get the information of the next vertex.
     *
     * @param vwgt The vertex weight(s).
     * @param list The adjacency list of the vertex.
     *
     * @return True if another vertex was found in the file.
     */
    virtual bool getNextVertex(
        std::vector<val_t> & vwgt,
        std::vector<MatrixEntry> & list) = 0;


    /**
     * @brief Get the name of the filetype.
     *
     * @return The name of the filetype.
     */
    virtual std::string const & getName() const = 0;


};




}




#endif
