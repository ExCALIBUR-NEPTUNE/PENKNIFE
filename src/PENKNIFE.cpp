///////////////////////////////////////////////////////////////////////////////
//
// File: main.cpp
//
//
// Description: Entrypoint for the PENKNIFE solver
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <mpi.h>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverCore/Driver/Driver.h>
#include <SpatialDomains/MeshGraphIO.h>

int main(int argc, char *argv[])
{

    // MPI is initialised/finalised here to ensure that Nektar++ does not
    // initialise/finalise MPI when we were not expecting it to.
    int provided_thread_level;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED,
                        &provided_thread_level) != MPI_SUCCESS)
    {
        std::cout << "ERROR: MPI_Init != MPI_SUCCESS" << std::endl;
        return -1;
    }
    int result = 0;
    try
    {
        auto session =
            Nektar::LibUtilities::SessionReader::CreateInstance(argc, argv);
        auto graph = Nektar::SpatialDomains::MeshGraphIO::Read(session);
        std::string module;
        session->LoadSolverInfo("Driver", module, "Standard");
        auto driver = Nektar::SolverCore::GetDriverFactory().CreateInstance(
            module, session, graph);
        driver->Execute();
        session->Finalise();
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        result = 1;
    }
    if (MPI_Finalize() != MPI_SUCCESS)
    {
        std::cout << "ERROR: MPI_Finalize != MPI_SUCCESS" << std::endl;
        return -1;
    }

    return result;
}
