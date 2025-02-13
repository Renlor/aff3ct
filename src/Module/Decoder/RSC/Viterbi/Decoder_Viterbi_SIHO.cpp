
#include "Module/Decoder/RSC/Viterbi/Decoder_Viterbi_SIHO.hpp"

#define DOUBLE_INF std::numeric_limits<double>::infinity()

using namespace aff3ct;
using namespace aff3ct::module;

template<typename B, typename R>
Decoder_Viterbi_SIHO<B, R>::Decoder_Viterbi_SIHO(const int K,
                                                 const std::vector<std::vector<int>> trellis,
                                                 const bool is_closed)
  : Decoder_SIHO<B, R>(K, is_closed ? 2 * K + 2 * static_cast<int>(std::log2(trellis[0].size())) : 2 * K)
  , m_n_states(static_cast<int>(trellis[0].size()))
  , m_n_memories(static_cast<int>(std::log2(m_n_states)))
  , m_K(K)
  , m_N(is_closed ? 2 * K + 2 * m_n_memories : 2 * K)
  , m_n_steps(K + m_n_memories)
  , m_trellis(trellis)
  , m_T({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_T_inv({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_C({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_P(std::vector<double>(m_n_states * (m_n_steps + 1)))
  , m_backwards_path(std::vector<int>(m_n_states * (m_n_steps + 1)))
  , m_M(std::vector<int>(m_n_states * m_n_states))
  , m_decoded(std::vector<int>(m_K + m_n_memories))
  , m_bin_vals(std::vector<std::array<double, 2>>{ { 0.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 }, { 1.0, 0.0 } })
  , m_closing_path(std::vector<int>(m_n_states))
  , m_closing_inputs(std::vector<int>(m_n_states))
{
    const std::string name = "Decoder_Viterbi_SIHO";
    this->set_name(name);
    for (auto& t : this->tasks)
        t->set_replicability(true);

    this->setup();
}

template<typename B, typename R>
Decoder_Viterbi_SIHO<B, R>*
Decoder_Viterbi_SIHO<B, R>::clone() const
{
    auto m = new Decoder_Viterbi_SIHO(*this);
    m->deep_copy(*this);
    return m;
}

/**
 * @brief Create all the matrices used by the Viterbi algorithm.
 *
 * @tparam B Output integer data type
 * @tparam R Input floating-point data type
 */
template<typename B, typename R>
void
Decoder_Viterbi_SIHO<B, R>::setup()
{
    // Filling T and C
    // T is the transition matrix, used as follow:
    //   next_state = m_T[sys_bit][previous_state]
    // C is the output matrix. It is used the same as T
    // but gives the output bits of the encoder for a given
    // transition.
    std::copy(m_trellis[6].begin(), m_trellis[6].end(), m_T[0].begin());
    std::copy(m_trellis[8].begin(), m_trellis[8].end(), m_T[1].begin());

    std::copy(m_trellis[7].begin(), m_trellis[7].end(), m_C[0].begin());
    std::copy(m_trellis[9].begin(), m_trellis[9].end(), m_C[1].begin());

    // Filling M, the input bit matrix.
    // For each possible transition of the encoder, M stores the associated
    // input bit + 1. An impossible transition is stored as a 0.
    int next_state;
    for (auto previous_state = 0; previous_state < m_n_states; previous_state++)
    {
        for (auto input_bit = 0; input_bit < 2; input_bit++)
        {
            next_state = m_T[input_bit][previous_state];
            m_M[previous_state * m_n_states + next_state] = input_bit + 1;
        }
    }

    // Inversion of T
    // This inversion gives the previous state given the next state
    for (auto bit_sys = 0; bit_sys < 2; bit_sys++)
    {
        for (auto state = 0; state < m_n_states; state++)
        {
            m_T_inv[bit_sys][m_T[bit_sys][state]] = state;
        }
    }

    // Filling closing paths and their associated inputs
    // Used for trellis termination
    for (auto i_state = 0; i_state < m_n_states; i_state++)
        if (m_T[0][i_state] < m_T[1][i_state])
        {
            m_closing_path[i_state] = m_T[0][i_state];
            m_closing_inputs[i_state] = 0;
        }
        else
        {
            m_closing_path[i_state] = m_T[1][i_state];
            m_closing_inputs[i_state] = 1;
        }
}

/**
 * @brief Reset the trellis before decoding
 *
 * @tparam B Output integer data type
 * @tparam R Input floating-point data type
 */
template<typename B, typename R>
void
Decoder_Viterbi_SIHO<B, R>::__reset()
{
    std::fill(m_P.begin(), m_P.end(), DOUBLE_INF);
    std::fill(m_backwards_path.begin(), m_backwards_path.end(), 0);
    std::fill(m_decoded.begin(), m_decoded.end(), 0);
}

// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
/**
 * @brief Test the sign of a value
 *
 * @tparam T type of the value
 * @param val value to test
 * @return int 1  if the value is positive
 *             0  if it is 0
 *             -1 if the value is negative
 */
template<typename T>
int
sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief Forward pass in the trellis. Calculate all branch and node metrics and create paths
 * in the trellis
 *
 * @tparam B Output integer data type
 * @tparam R Input floating-point data type
 * @param Y_N Input LLR pointer
 */
template<typename B, typename R>
void
Decoder_Viterbi_SIHO<B, R>::_forward_pass(const R* Y_N)
{
    // Initialize first state at node metric 0
    m_P[0] = 0.0;

    for (auto i_step = 1; i_step <= m_K; i_step++)
    {
        // LLR extraction
        std::array<R, 2> channel_input{ Y_N[2 * (i_step - 1)], Y_N[2 * (i_step - 1) + 1] };

        // LLR saturation
        for (auto i = 0; i < 2; i++)
        {
            if (channel_input[i] > 100) channel_input[i] = 100;
            if (channel_input[i] < -100) channel_input[i] = -100;
        }

        for (auto i_state = 0; i_state < m_n_states; i_state++)
        {
            for (auto bit_sys = 0; bit_sys < 2; bit_sys++)
            {
                const int previous_state = m_T_inv[bit_sys][i_state];
                const int previous_node_idx = (i_step - 1) * m_n_states + previous_state;
                const int next_node_idx = (i_step)*m_n_states + i_state;
                // Si on peut partir de l'etat donne a l'etape donnee
                if (m_P[previous_node_idx] != DOUBLE_INF)
                {
                    // Recuperation du poids du noeud
                    const double node_weight = m_P[previous_node_idx];

                    const std::array<double, 2> output = m_bin_vals[2 * bit_sys + m_C[bit_sys][previous_state]];

                    const double new_weight = node_weight + channel_input[0] * output[0] + channel_input[1] * output[1];

                    // s'il existe une etape suivante
                    if (new_weight < m_P[next_node_idx])
                    {
                        m_P[next_node_idx] = new_weight;
                        m_backwards_path[next_node_idx] = previous_state;
                    }
                }
            }
        }
    }
    for (auto i_step = m_K; i_step < m_K + m_n_memories; i_step++)
    {
        std::array<R, 2> channel_input{ Y_N[2 * i_step], Y_N[2 * i_step + 1] };

        for (auto i = 0; i < 2; i++)
        {
            if (channel_input[i] > 100) channel_input[i] = 100;
            if (channel_input[i] < -100) channel_input[i] = -100;
        }

        for (auto i_state = 0; i_state < m_n_states; i_state++)
        {
            const int node_idx = i_step * m_n_states + i_state;

            if (m_P[node_idx] != DOUBLE_INF)
            {
                const double node_weight = m_P[node_idx];
                const int next_state = m_closing_path[i_state];
                const int bit_sys = m_closing_inputs[i_state];
                const int next_node_idx = (i_step + 1) * m_n_states + next_state;

                const std::array<double, 2> output = m_bin_vals[2 * bit_sys + m_C[bit_sys][i_state]];

                const double new_weight = node_weight + channel_input[0] * output[0] + channel_input[1] * output[1];

                if (new_weight < m_P[next_node_idx])
                {
                    m_P[next_node_idx] = new_weight;
                    m_backwards_path[next_node_idx] = i_state;
                }
            }
        }
    }
}

/**
 * @brief Backwards pass for the Viterbi algorithm.
 * Follow the best path along the trellis to decode the received message.
 * Safely assumes the decoding starts at state 0
 *
 * @tparam B
 * @tparam R
 * @param V_K Output pointer for the decoded message
 */
template<typename B, typename R>
void
Decoder_Viterbi_SIHO<B, R>::_backwards_pass(B* V_K)
{
    int state = 0;
    for (auto i_step = m_K + m_n_memories - 1; i_step >= 0; i_step--)
    {
        const int next_state = state;
        const int previous_state = m_backwards_path[(i_step + 1) * m_n_states + next_state];

        if (i_step < m_K)
        {
            V_K[i_step] = m_M[previous_state * m_n_states + next_state] - 1;
        }
        state = previous_state;
    }
}

template<typename B, typename R>
int
Decoder_Viterbi_SIHO<B, R>::_decode_siho(const R* Y_N, B* V_K, const size_t /*frame_id*/)
{
    __reset();
    _forward_pass(Y_N);
    _backwards_pass(V_K);

    return 0;
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Decoder_Viterbi_SIHO<B_8, Q_8>;
template class aff3ct::module::Decoder_Viterbi_SIHO<B_16, Q_16>;
template class aff3ct::module::Decoder_Viterbi_SIHO<B_32, Q_32>;
template class aff3ct::module::Decoder_Viterbi_SIHO<B_64, Q_64>;
#else
template class aff3ct::module::Decoder_Viterbi_SIHO<B, Q>;
#endif
// ==================================================================================== explicit template instantiation
