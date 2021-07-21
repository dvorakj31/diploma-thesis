#include "Lenstra.h"
#include <iostream>
#include <boost/mpi.hpp>
#include <omp.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace mpi = boost::mpi;

NTL::ZZ generated_counter{0};
int end = 0;

NTL::ZZ Lenstra::factorize() const {
    /// Sequential algorithm for computing factorization
    auto divisor = NTL::conv<NTL::ZZ>(1);
    auto sqrt_n = NTL::SqrRoot(*_options->composite_number);

    NTL::ZZ bound = sqrt_n;
    if (*_options->bound > 2)
        bound = (*_options->bound < sqrt_n ? *_options->bound : sqrt_n);
    auto divided = bound / 1000000;
    NTL::ZZ test_after{(divided < 100 ? (NTL::ZZ)100 : divided)};
    NTL::ZZ counter(1);
    while (true) {
        auto point = _model->generate_elliptic_curve();
        for (NTL::ZZ k{2}; k < bound; k++, counter++) {
            point = _model->mul_points(k, point);
            if (counter % test_after == 0) {
                divisor = _model->try_get_factor(point);
                if (divisor > 1 && divisor < *_options->composite_number) {
                    return divisor;
                }
                counter = 0;
            }
            if (_model->is_infinity_point(point)) {
                break;
            }
        }
        divisor = _model->try_get_factor(point);
        if (divisor > 1 && divisor < *_options->composite_number) {
            return divisor;
        }
    }
}

NTL::ZZ Lenstra::factorize_parallel(int argc, char **argv) {
    /// Master-slave algorithm
    /// First we initialize communication environment via mpi::environment and mpi::communicator
    mpi::environment env(argc, argv, mpi::threading::multiple);
    mpi::communicator world;

    /// Result of factorizing
    NTL::ZZ result{0};

    if (!world.rank()) {
        /// master part
        _generate_ecc(env, world);
        _point = _model->generate_elliptic_curve();
    } else {
        /// slave part
        _get_ecc(env, world);
    }

    mpi::timer timer;

    result = _factorize_parallel(env, world);

    double elapsed = timer.elapsed();

    if (result != 0) {
        std::cout << "proc " << world.rank() << ": factor = " << result << "\n";
        std::cout << "proc " << world.rank() << ": time = " << elapsed << " s\n";
        _stop_all(env, world);
        mpi::environment::abort(0);
    }

    world.barrier();

    if (world.rank() == 0) {
        std::cout << "generated ecc = " << generated_counter << "\n";
    }

    return result;
}

void Lenstra::_generate_ecc(const mpi::environment &environment, const mpi::communicator &communicator) {
    /// Generates elliptic curve for all working processes
    for (int i = 1; i < communicator.size(); i++) {
        _par_generate_ecc(environment, communicator);
    }
}

NTL::ZZ Lenstra::_factorize_parallel(const mpi::environment &environment, const mpi::communicator &communicator) {
    /// starts computing algorithm for process
    std::vector<ProjectivePoint> points;
    NTL::ZZ result{0};
    NTL::ZZ k{1};
    NTL::ZZ counter{1};
    auto divisor = NTL::conv<NTL::ZZ>(1);
    auto sqrt_n = NTL::SqrRoot(*_options->composite_number);

    NTL::ZZ bound = sqrt_n;
    if (*_options->bound > 2)
        bound = (*_options->bound < sqrt_n ? *_options->bound : sqrt_n);
    auto divided = bound / 1000000;
    NTL::ZZ test_after{(divided < 100 ? (NTL::ZZ)100 : divided)};

    while(!end) {
        k = 2;
        points.clear();
        #pragma omp parallel shared(points, result, end, k, environment, communicator, sqrt_n, test_after)
        {
            NTL::ZZ tmp;
            auto point = _point;
            bool tested = false;
            while (k < bound && !end) {
                #pragma omp critical
                {
                    tmp = k;
                    k++;
                    counter++;
                }

                point = _model->mul_points(tmp, point);
                if (_model->is_infinity_point(point)) {
                    break;
                }

                if (!tested && counter >= test_after) {
                    tested = true;
                    auto factor = _model->try_get_factor(point);
                    if (factor > 1 && factor < *_options->composite_number) {
                        result = factor;
                        #pragma omp critical (ending)
                        {
                            end = 1;
                        }
                        break;
                    }
                }

                #pragma omp critical (ending)
                {
                    if (!end && communicator.rank() == 0) {
                        end = !_par_generate_ecc(environment, communicator);
                    } else if (!end && communicator.rank() != 0) {
                        end = _check_end(environment, communicator) || end;
                    }
                }
            }

            #pragma omp critical (add_point)
            {
                points.push_back(point);
            }

            #pragma omp barrier

            if (counter % test_after == 0) {
                auto factor = _model->try_get_factor(point);
                if (factor > 1 && factor < *_options->composite_number) {
                    result = factor;
                    end = 1;
                }
            }

            #pragma omp single
            {
                #pragma omp critical (ending)
                {
                    if (communicator.rank() != 0)
                        end = _check_end(environment, communicator);
                    else
                        end = !_par_generate_ecc(environment, communicator);
                }
                for (const auto &vec_point : points) {
                    _point = _model->add_points(_point, vec_point);
                }

                if (!_model->is_infinity_point(_point)) {
                    auto factor = _model->try_get_factor(_point);
                    if (factor > 1 && factor < *_options->composite_number) {
                        result = factor;
                        end = 1;
                    }
                }
                if (!end) {
                    if (communicator.rank() == 0) {
                        generated_counter++;
                        _point = _model->generate_elliptic_curve();
                        end = !_par_generate_ecc(environment, communicator);
                    } else {
                        end = !_get_ecc(environment, communicator);
                    }
                }
            }
        }
    }
    return result;
}

void Lenstra::_generate_edwards(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator,
                           int source) {
    /// Generates edwards curve  for working process
    auto model = dynamic_cast<EdwardsModel*>(_model.get());
    std::stringstream buffer;
    boost::archive::text_oarchive ar(buffer);
    auto point = model->generate_elliptic_curve();
    auto curve = model->get_elliptic_curve();

    ar << point << curve;
    generated_counter++;
    communicator.send(source, TAGS::NEW_ECC, buffer.str());
}

void Lenstra::_generate_weierstrass(const boost::mpi::environment &environment, const boost::mpi::communicator &communicator,
                               int source) {
    /// Generates weierstrass curve for working process
    auto model = dynamic_cast<WeierstrassModel*>(_model.get());
    std::stringstream buffer;
    boost::archive::text_oarchive ar(buffer);
    auto point = model->generate_elliptic_curve();
    auto curve = model->get_elliptic_curve();

    ar << point << curve;
    generated_counter++;
    communicator.send(source, TAGS::NEW_ECC, buffer.str());
}

bool Lenstra::_get_ecc(const mpi::environment &environment, const mpi::communicator &communicator) {
    /// Gets new elliptic curve from master process
    communicator.send(0, TAGS::NEW_ECC);
    auto status = communicator.probe();
    if (status.tag() == TAGS::NEW_ECC) {
        std::string msg;
        communicator.recv(status.source(), status.tag(), msg);
        std::istringstream buffer(msg);
        boost::archive::text_iarchive ar(buffer);

        ar >> _point;
        _point.z = 1;
        if (_options->weierstrass) {
            WeierstrassModel::EllipticCurve ecc;
            ar >> ecc;
            dynamic_cast<WeierstrassModel*>(_model.get())->set_elliptic_curve(ecc);

        } else {
            EdwardsModel::EllipticCurve ecc;
            ar >> ecc;
            dynamic_cast<EdwardsModel*>(_model.get())->set_elliptic_curve(ecc);
        }
        return true;
    }
    return false;
}

bool Lenstra::_check_end(const mpi::environment &environment, const mpi::communicator &communicator) {
    /// Check if we have pending message
    auto status = communicator.iprobe();
    if (status.has_value()) {
        communicator.recv(status.value().source(), status.value().tag());
    }
    return status.has_value() && status.value().tag() == TAGS::STOP;
}

void Lenstra::_stop_all(const mpi::environment &, const boost::mpi::communicator& communicator) {
    /// Send stop tag to all processses
    std::vector<mpi::request> requests;
    for (int i = 0; i < communicator.size(); i++) {
        if (i != communicator.rank()) {
            requests.emplace_back(communicator.isend(i, TAGS::STOP));
        }
    }
    mpi::wait_all(requests.begin(), requests.end());
}

bool Lenstra::_par_generate_ecc(const mpi::environment &environment, const mpi::communicator &communicator) {
    /// Master part with generating of new elliptic curve
    auto status = communicator.iprobe();
    if (status.has_value() && status.value().tag() == TAGS::STOP) {
        communicator.recv(status.value().source(), status.value().tag());
        return false;
    } else if (status.has_value() && status.value().tag() == TAGS::NEW_ECC) {
        communicator.recv(status.value().source(), status.value().tag());
        if (_options->weierstrass) {
            _generate_weierstrass(environment, communicator, status.value().source());
        } else {
            _generate_edwards(environment, communicator, status.value().source());
        }
    }
    return true;
}
