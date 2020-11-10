/*
 * Abstraction.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/** @file **/
#ifndef ABSTRACTION_HH_
#define ABSTRACTION_HH_

#include <iostream>
#include <fstream>
#include <cstring>
#include <memory>
#include <vector>

#include "UniformGrid.hh"
#include "TransitionFunction.hh"
#include "json.hpp"

/** @namespace scots **/
namespace scots {

/** @cond **/
/* default parameter for the third parameter of Abstraction::compute_gb  */
    namespace params {
        auto avoid_abs = [](const abs_type&) noexcept {return false;};
    }
/** @endcond **/


/**
 * @class Abstraction
 * 
 * @brief Computation of the transition function of symbolic model/abstraction
 *
 * It implements in compute_gb the computation of the transition function based
 * on a growth bound. Additionally, it provides the functions print_post and
 * print_post_gb to print the center of the cells that cover the attainable set
 * which is computed from the growth bound.
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual.pdf</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 *
 **/
    template<class state_type, class input_type>
    class Abstraction {
    private:
        /* grid information of state alphabet */
        const UniformGrid m_state_alphabet;
        /* grid information of input alphabet */
        const UniformGrid m_input_alphabet;
        /* measurement error bound */
        std::unique_ptr<double[]> m_z;
        /* print progress to the console (default m_verbose=true) */
        bool m_verbose=true;
        /* to display the progress of the computation of the abstraction */
        void progress(const abs_type& i, const abs_type& N, abs_type& counter) {
            if(!m_verbose)
                return;
            if(((double)i/(double)N*100)>=counter){
                if((counter%10)==0)
                    std::cout << counter;
                else if((counter%2)==0)
                    std::cout << ".";
                counter++;
            }
            std::flush(std::cout);
            if(i==(N-1))
                std::cout << "100\n";
        }
    public:
        /* json path for reading results from cora */
        std::string cora_path;
        /* json path for writing comparison results */
        std::string compare_path;
        /* @cond  EXCLUDE from doxygen*/
        /* deactivate standard constructor */
        Abstraction() = delete;
        /* cannot be copied or moved */
        Abstraction(Abstraction&&) = delete;
        Abstraction(const Abstraction&) = delete;
        Abstraction& operator=(Abstraction&&)=delete;
        Abstraction& operator=(const Abstraction&)=delete;
        /* @endcond */

        /** @brief constructor with the abstract state alphabet and abstract input
         * alphabet (measurement error bound is set to zero)
         *
         *  @param state_alphabet - UniformGrid representing the abstract state alphabet
         *  @param input_alphabet - UniformGrid representing the abstract input alphabet
         **/
        Abstraction(const UniformGrid& state_alphabet,
                    const UniformGrid& input_alphabet) :
                m_state_alphabet(state_alphabet),
                m_input_alphabet(input_alphabet),
                m_z(new double[state_alphabet.get_dim()]()) {
            /* default value of the measurement error
               * (heurisitc to prevent rounding errors)*/
            for(int i=0; i<m_state_alphabet.get_dim(); i++)
                m_z[i]=m_state_alphabet.get_eta()[i]/1e10;
        }

        /**
         * @brief computes the transition function using linearization
         *
         * @param[out] transition_function - the result of the computation
         *
         * @param[in] system_post - lambda expression of the form
         *                          \verbatim [] (state_type &x, const input_type &u) ->  void  \endverbatim
         *                          system_post(x,u) provides a numerical approximation of ODE
         *                          solution at time tau with initial state x and input u \n
         *                          the result is stored in x
         *
         * @param[in] radius_post - lambda expression of the form
         *                          \verbatim [] (state_type &r, const state_type& x, const input_type &u) -> void  \endverbatim
         *                          radius_post(x,u) provides a numerical approximation of
         *                          the growth bound for the cell (with center x, radius  r) and input u\n
         *                          the result is stored in r
         *
         * @param[in] avoid  - OPTIONALLY provide lambda expression of the form
         *                     \verbatim [] (const abs_type &i) -> bool \endverbatim
         *                     returns true if the abstract state i is in the avoid
         *                     set; otherwise returns false
         *
         * @param[in] gb - OPTIONALLY flag is true if the reachable set is computed using the growth bound approach.
         * Otherwise, data will be read from CORA json in the @cora_path.
         *
         * @param[in] compare - OPTIONALLY
         *
         * The computation proceeds in two loops. In the first loop the cornerIDs are
         * computed, which represent the cell IDs that cover the over-approximation of the attainable set:
         * \verbatim corner_IDs[i*M+j+0] \endverbatim = lower-left cell ID of the
         * integer hyper-interval of cell IDs that cover the over-approximation of the
         * attainable set of cell with ID=i under input ID=j
         * \verbatim corner_IDs[i*M+j+1] \endverbatim = upper-right  of the
         * integer hyper-interval of cell IDs that cover the over-approximation of the
         * attainable set of cell with ID=i under input ID=j
         *
         * In the second loop the data members of the TransitionFunction are computed.
         *
         **/
        template<class F1, class F2, class F3=decltype(params::avoid_abs)>
        void compute_gb(TransitionFunction& transition_function,
                        F1& system_post,
                        F2& radius_post,
                        F3& avoid=params::avoid_abs,
                        bool gb=true,
                        bool compare = false) {
            /* number of cells */
            abs_type N=m_state_alphabet.size();
            /* number of inputs */
            abs_type M=m_input_alphabet.size();
            /* number of transitions (to be computed) */
            abs_ptr_type T=0;
            /* state space dimension */
            int dim=m_state_alphabet.get_dim();
            /* for display purpose */
            abs_type counter=0;
            /* some grid information */
            std::vector<abs_type> NN=m_state_alphabet.get_nn();
            /* variables for managing the post */
            std::vector<abs_type> lb(dim);  /* lower-left corner */
            std::vector<abs_type> ub(dim);  /* upper-right corner */
            std::vector<abs_type> no(dim);  /* number of cells per dim */
            std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
            /* radius of hyper interval containing the attainable set */
            state_type eta;
            state_type r;
            /* state and input variables */
            state_type x;
            input_type u;
            /* state and input variables for logging results in json */
            state_type state_temp;
            input_type input_temp;
            /* for out of bounds check */
            state_type lower_left;
            state_type upper_right;
            /* for storing comparison results with cora */
            nlohmann::json json_comparison = nlohmann::json::array();
            std::ofstream myfile;
            /* for reading data from cora */
            nlohmann::json::iterator it;
            nlohmann::json json_cora;
            /* initialization for reading data from cora or comparison */
            if((compare && gb) || !gb){
                std::ifstream ipstr(cora_path);
                ipstr >> json_cora;
                it = json_cora.begin();
                if(compare){
                    myfile.open(compare_path);
                }
            }

            /* copy data from m_state_alphabet */
            for(int i=0; i<dim; i++) {
                eta[i]=m_state_alphabet.get_eta()[i];
                lower_left[i]=m_state_alphabet.get_lower_left()[i];
                upper_right[i]=m_state_alphabet.get_upper_right()[i];
            }
            /* init in transition_function the members no_pre, no_post, pre_ptr */
            transition_function.init_infrastructure(N,M);
            /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
            std::unique_ptr<abs_type[]> corner_IDs(new abs_type[N*M*2]());
            /* is post of (i,j) out of domain ? */
            std::unique_ptr<bool[]> out_of_domain(new bool[N*M]());
            /*
             * first loop: compute corner_IDs:
             * corner_IDs[i*M+j][0] = lower-left cell index of over-approximation of attainable set
             * corner_IDs[i*M+j][1] = upper-right cell index of over-approximation of attainable set
             */
            /* loop over all states. TODO: Try to parallelize this loop */
            for(abs_type i=0; i<N; i++) {
                /* is i an element of the avoid symbols ? */
                if(avoid(i)) {
                    for(abs_type j=0; j<M; j++) {
                        out_of_domain[i*M+j]=true;
                    }
                    continue;
                }
                /* loop over all inputs */
                for(abs_type j=0; j<M; j++) {
                    out_of_domain[i*M+j]=false;
                    /* get center x of cell */
                    m_state_alphabet.itox(i,x);
                    /* cell radius (including measurement errors) */
                    for(int k=0; k<dim; k++) // TODO: put it on the top
                        r[k]=eta[k]/2.0+m_z[k];
                    /* current input */
                    m_input_alphabet.itox(j,u);

                    /* determine the cells which intersect with the attainable set:
                     * discrete hyper interval of cell indices
                     * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
                     * covers attainable set
                     */
                    abs_type npost=1;

                    /* integrate system and radius growth bound for provided input and state. */
                    /* the result is stored in x and r */
                    if(gb) {
                        if(compare){
                            state_temp = x;
                            input_temp = u;
                        }
                        radius_post(r, x, u); // solution with disturbance
                        system_post(x, u); // nominal solution
                    }

                    /* initialization for jsons */
                    auto left_scots = nlohmann::json::array();
                    auto right_scots = nlohmann::json::array();
                    auto left_cora = nlohmann::json::array();
                    auto right_cora = nlohmann::json::array();
                    auto no_cora = nlohmann::json::array();
                    auto no_scots = nlohmann::json::array();
                    /* left and right bounds of interval for each dimension */
                    double left;
                    double right;
                    double lb_cora_k;
                    double ub_cora_k;
                    /* find bounds for all dimensions. */
                    for (int k = 0; k < dim; k++) { // TODO: re-check that dimension's loop is right
                        if (gb) {
                            left = x[k] - r[k] - m_z[k];
                            right = x[k] + r[k] + m_z[k];
                        }
                        if ((compare && gb) || !gb) {
                            nlohmann::json j1 = *it;
                            auto cen = j1.at("zonotope").at("value").at("center")[k].get<double>();
                            auto g = j1.at("zonotope").at("value")["generators"][k][k].get<double>();
                            ++it;

                            if (compare) {
                                lb_cora_k = cen - g - m_z[k];
                                ub_cora_k = cen + g - m_z[k];
                                left_cora.push_back(lb_cora_k);
                                right_cora.push_back(ub_cora_k);
                                left_scots.push_back(left);
                                right_scots.push_back(right);
                            } else {
                                left = cen - g - m_z[k];
                                right = cen + g + m_z[k];
                            }
                        }

                        /* check for out of bounds */
                        if (left <= lower_left[k] - eta[k] / 2.0 || right >= upper_right[k] + eta[k] / 2.0) {
                            out_of_domain[i * M + j] = true;
                            std::cout << "out of domain " << '\n';
                            break;
                        }

                        /* update bounds. number of steps
                         * integer coordinate of lower left corner of post
                         * */
                        lb[k] = static_cast<abs_type>((left - lower_left[k] + eta[k] / 2.0) / eta[k]);
                        /* integer coordinate of upper right corner of post */
                        ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
                        /* number of grid points in the post in each dimension */
                        no[k] = (ub[k] - lb[k] + 1);
                        if(compare){
                            no_scots.push_back(no[k]);

                            no_cora.push_back((lb_cora_k- lower_left[k] + eta[k]/2.0)/eta[k]-
                                              (ub_cora_k-lower_left[k]+eta[k]/2.0)/eta[k] + 1);
                        }
                        /* total number of post */
                        npost *= no[k];
                        cc[k] = 0;
                    }
                    if (compare) {
                        nlohmann::json j;
                        j["zonotope"] = {"input", input_temp, "state", state_temp,
                                         "scots", {"lb", left_scots, "ub", right_scots, "steps", no_scots},
                                         "cora", {"lb", left_cora, "ub", right_cora, "steps", no_cora}};

                        json_comparison.push_back(j);
                    }

                    corner_IDs[i*(2*M)+2*j]=0;
                    corner_IDs[i*(2*M)+2*j+1]=0;
                    if(out_of_domain[i*M+j])
                        continue;

                    /* compute indices of post */
                    for(abs_type k=0; k<npost; k++) {
                        abs_type q=0;
                        for(int l=0; l<dim; l++)
                            q+=(lb[l]+cc[l])*NN[l];
                        cc[0]++;
                        for(int l=0; l<dim-1; l++) {
                            if(cc[l]==no[l]) {
                                cc[l]=0;
                                cc[l+1]++;
                            }
                        }
                        /* (i,j,q) is a transition */
                        /* increment number of pres for (q,j) */
                        transition_function.m_no_pre[q*M+j]++;
                        /* store id's of lower-left and upper-right cell */
                        if(k==0)
                            corner_IDs[i*(2*M)+2*j]=q;
                        if(k==npost-1)
                            corner_IDs[i*(2*M)+2*j+1]=q;
                    }
                    /* increment number of transitions by number of post */
                    T+=npost;
                    transition_function.m_no_post[i*M+j]=npost;
                }
                /* print progress */
                if(m_verbose) {
                    if(counter==0)
                        std::cout << "1st loop: ";
                }
                progress(i,N,counter);
            }

            /* write comparison result with cora to the file */
            if (compare) {
                std::ofstream myfile;
                myfile.open(compare_path);
                myfile << json_comparison;
                myfile.close();
            }

            /* compute pre_ptr */
            abs_ptr_type sum=0;
            for(abs_type i=0; i<N; i++) {
                for(abs_type j=0; j<M; j++) {
                    sum+=transition_function.m_no_pre[i*M+j];
                    transition_function.m_pre_ptr[i*M+j]=sum;
                }
            }
            /* allocate memory for pre list */
            transition_function.init_transitions(T);

            /* second loop: fill pre array */
            counter=0;
            for(abs_type i=0; i<N; i++) {
                /* loop over all inputs */
                for(abs_type j=0; j<M; j++) {
                    /* is x an element of the overflow symbols ? */
                    if(out_of_domain[i*M+j])
                        continue;
                    /* extract lower-left and upper-bound points */
                    abs_type k_lb=corner_IDs[i*2*M+2*j];
                    abs_type k_ub=corner_IDs[i*2*M+2*j+1];
                    abs_type npost=1;

                    /* cell idx to coordinates */
                    for(int k=dim-1; k>=0; k--) {
                        /* integer coordinate of lower left corner */
                        lb[k]=k_lb/NN[k];
                        k_lb=k_lb-lb[k]*NN[k];
                        /* integer coordinate of upper right corner */
                        ub[k]=k_ub/NN[k];
                        k_ub=k_ub-ub[k]*NN[k];
                        /* number of grid points in each dimension in the post */
                        no[k]=(ub[k]-lb[k]+1);
                        /* total no of post of (i,j) */
                        npost*=no[k];
                        cc[k]=0;

                    }

                    for(abs_type k=0; k<npost; k++) {
                        abs_type q=0;
                        for(int l=0; l<dim; l++)
                            q+=(lb[l]+cc[l])*NN[l];
                        cc[0]++;
                        for(int l=0; l<dim-1; l++) {
                            if(cc[l]==no[l]) {
                                cc[l]=0;
                                cc[l+1]++;
                            }
                        }
                        /* (i,j,q) is a transition */
                        transition_function.m_pre[--transition_function.m_pre_ptr[q*M+j]]=i;
                    }
                }
                /* print progress */
                if(m_verbose) {
                    if(counter==0)
                        std::cout << "2nd loop: ";
                }
                progress(i,N,counter);
            }
        }

        /**
         * @brief computes the transition function
         *
         * @param[out] transition_function - the result of the computation
         * @param[in] path -   provide path to read reachable set that computed using linearization technique
         *                     from https://ieeexplore.ieee.org/document/6632887. Otherwise, the growth bound approach
         *                     will be used.
         *
         * @param[in] system_post - same as in @compute_gb.
         *
         * @param[in] radius_post - same as in @compute_gb.
         *
         * @param[in] avoid  - same as in @compute_gb.
         *
         * The computation proceeds in two loops as in growth bound approach. The only difference is in the first loop:
         * overapproximation of the reachable set is computed using linearization approach, not the growth bound.
         * Reachable set is computed using CORA framework - https://github.com/TUMcps/CORA.
         *
         **/
        template<class F1, class F2, class F3=decltype(params::avoid_abs)>
        void compute_transitions(TransitionFunction& transition_function,
                                   std::string path,
                                   F1& system_post,
                                   F2& radius_post,
                                   F3& avoid=params::avoid_abs){
            cora_path = path;
            compute_gb(transition_function, system_post, radius_post, params::avoid_abs, false, false);
        }

        /**
         * @brief computes the gb and compare results with cora
         *
         * @param[out] transition_function - the result of the computation
         *
         * @param[in] system_post - same as in @compute_gb.
         *
         * @param[in] radius_post - same as in @compute_gb.
         *
         * @param[in] avoid  - same as in @compute_gb.
         * @param[in] path -    provide path to read reachable set that computed using linearization technique
         *                     from https://ieeexplore.ieee.org/document/6632887. Otherwise, the growth bound approach
         *                     will be used.
         *
         * The computation proceeds in two loops as in growth bound approach. The only difference is in the first loop:
         * overapproximation of the reachable set is computed using linearization approach, not the growth bound.
         * Reachable set is computed using CORA framework - https://github.com/TUMcps/CORA.
         *
         **/
        template<class F1, class F2, class F3=decltype(params::avoid_abs)>
        void compute_compare_gb(TransitionFunction& transition_function,
                                 std::string path_cora,
                                 std::string path_compare,
                                 F1& system_post,
                                 F2& radius_post,
                                 F3& avoid=params::avoid_abs){
            cora_path = path_cora;
            compare_path = path_compare;
            compute_gb(transition_function, system_post, radius_post, params::avoid_abs, true, true);
        }

        /** @brief get the center of cells that are used to over-approximated the
         *  attainable set associated with cell (with center x) and input u
         *
         *  @returns std::vector<abs_type>
         *  @param system_post - lambda expression as in compute_gb
         *  @param radius_post - lambda expression as in compute_gb
         *  @param x - center of cell
         *  @param u - input
         **/
        template<class F1, class F2>
        std::vector<state_type> get_post(F1& system_post,
                                         F2& radius_post,
                                         const state_type& x,
                                         const input_type& u) const {
            /* initialize return value */
            std::vector<state_type> post {};
            /* state space dimension */
            int dim = m_state_alphabet.get_dim();
            /* variables for managing the post */
            std::vector<abs_type> lb(dim);  /* lower-left corner */
            std::vector<abs_type> ub(dim);  /* upper-right corner */
            std::vector<abs_type> no(dim);  /* number of cells per dim */
            std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
            /* get radius */
            state_type r;
            state_type eta;
            /* for out of bounds check */
            state_type lower_left;
            state_type upper_right;
            /* fill in data */
            for(int i=0; i<dim; i++) {
                eta[i]=m_state_alphabet.get_eta()[i];
                r[i]=m_state_alphabet.get_eta()[i]/2.0+m_z[i];
                lower_left[i]=m_state_alphabet.get_lower_left()[i];
                upper_right[i]=m_state_alphabet.get_upper_right()[i];
            }
            state_type xx=x;
            /* compute growth bound and numerical solution of ODE */
            // TODO: check where this is used
            radius_post(r,x,u);
            system_post(xx,u);

            /* determine post */
            abs_type npost=1;
            for(int k=0; k<dim; k++) {
                /* check for out of bounds */
                double left = xx[k]-r[k]-m_z[k];
                double right = xx[k]+r[k]+m_z[k];
                if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0) {
                    return post;
                }
                /* integer coordinate of lower left corner of post */
                lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
                /* integer coordinate of upper right corner of post */
                ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
                /* number of grid points in the post in each dimension */
                no[k]=(ub[k]-lb[k]+1);
                /* total number of post */
                npost*=no[k];
                cc[k]=0;
            }
            /* compute indices of post */
            for(abs_type k=0; k<npost; k++) {
                abs_type q=0;
                for(int l=0; l<dim; l++)  {
                    q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
                }
                cc[0]++;
                for(int l=0; l<dim-1; l++) {
                    if(cc[l]==no[l]) {
                        cc[l]=0;
                        cc[l+1]++;
                    }
                }
                /* (i,j,q) is a transition */
                m_state_alphabet.itox(q,xx);
                post.push_back(xx);
            }
            return post;
        }


        /** @brief print the center of cells that are stored in the transition
         *  function as post of the cell with center x and input u
         *
         *  @param transition_function - the transition function of the abstraction
         *  @param x - center of cell
         *  @param u - input
         **/
        void print_post(const TransitionFunction& transition_function,
                        const state_type& x,
                        const input_type& u) const {
            std::vector<state_type> post = get_post(transition_function, x, u);
            if(!post.size())
                std::cout << "\nPost is out of domain\n";
            std::cout << "\nPost states: \n";
            for(abs_type v=0; v<post.size(); v++) {
                for(int i=0; i<m_state_alphabet.get_dim(); i++) {
                    std::cout << post[v][i] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        /** @brief print the center of cells that are used to over-approximated the
         *  attainable set (computed according to the growth bound)
         *  associated with the cell with center x and input u
         *
         *  @param system_post - lambda expression as in compute_gb
         *  @param radius_post - lambda expression as in compute_gb
         *  @param x - center of cell
         *  @param u - input
         **/
        template<class F1, class F2>
        void print_post_gb(F1& system_post,
                           F2& radius_post,
                           const state_type& x,
                           const input_type& u) const {
            std::vector<state_type> post = get_post(system_post,radius_post,x,u);
            if(!post.size()) {
                std::cout << "\nPost is out of domain\n";
                return;
            }

            std::cout << "\nPost states: \n";
            for(abs_type v=0; v<post.size(); v++) {
                for(int i=0; i<m_state_alphabet.get_dim(); i++) {
                    std::cout << post[v][i] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        /** @brief set the measurement error bound **/
        void set_measurement_error_bound(const state_type& error_bound) {
            for(int i=0; i<m_state_alphabet.get_dim(); i++) {
                m_z[i]=error_bound[i];
            }
        }
        /** @brief get measurement error bound **/
        std::vector<double> get_measruement_error_bound() {
            std::vector<double> z;
            for(int i=0; i<m_state_alphabet.get_dim(); i++) {
                z.push_back(m_z[i]);
            }
            return z;
        }
        /** @brief activate console output **/
        void verbose_on() {
            m_verbose=true;
        }
        /** @brief deactivate console output **/
        void verbose_off() {
            m_verbose=false;
        }
    };

} /* close namespace */
#endif /* ABSTRACTION_HH_ */
