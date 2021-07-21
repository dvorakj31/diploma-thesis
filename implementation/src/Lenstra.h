#ifndef DIP_LENSTRA_H
#define DIP_LENSTRA_H

#include <NTL/ZZ.h>

#include <utility>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "AbstractModel.h"
#include "EdwardsModel.h"
#include "WeierstrassModel.h"

class Lenstra final {
public:
    explicit Lenstra(std::shared_ptr<Options> options, std::shared_ptr<AbstractModel> model)
        : _options(std::move(options)), _model(std::move(model)) {
    }

    /// This function is used for sequential computation of factor
    [[nodiscard]] NTL::ZZ factorize() const;

    /// This function is used for parallel computation of factor
    [[nodiscard]] NTL::ZZ factorize_parallel(int argc, char **argv);

private:

    /// TAGS for signalizing processes what to do.
    enum TAGS {
        NEW_ECC = 0x1000,
        STOP = 0x0100,
    };

    std::shared_ptr<Options> _options;
    std::shared_ptr<AbstractModel> _model;
    /// Stores point for parallel purpose
    ProjectivePoint _point;

    /// This method is for process computation. It uses OpenMP pragmas.
    NTL::ZZ _factorize_parallel(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator);
    /// Auxiliary function for master process. Generates new elliptic curve.
    void _generate_ecc(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator);
    /// Auxiliary function for master process. Generates and sends new Edwards curve to working process specified by source argument.
    void _generate_edwards(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator,
                           int source);
    /// Auxiliary function for master process. Generates and sends new Weierstrass curve to working process specified by source argument.
    void _generate_weierstrass(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator,
                               int source);
    /// Auxiliary function for working process. Gets new elliptic curve for working process.
    bool _get_ecc(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator);

    /// Auxiliary function for checking if some process finished it's job.
    bool _check_end(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator);

    /// This function sends message to all other processes with end indication
    void _stop_all(const boost::mpi::environment &, const boost::mpi::communicator& communicator);

    /// Auxiliary function for
    bool _par_generate_ecc(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator);
};


#endif //DIP_LENSTRA_H
