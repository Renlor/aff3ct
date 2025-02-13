#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "Module/Decoder/RSC/Viterbi_list/Decoder_Viterbi_list_parallel.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

const double gDOUBLE_INF = std::numeric_limits<double>::infinity();

template<typename B, typename Q>
Decoder_Viterbi_list_parallel<B, Q>::Decoder_Viterbi_list_parallel(const int K,
                                                                   const int N,
                                                                   const int L,
                                                                   const module::CRC<B>& crc,
                                                                   const std::vector<std::vector<int>> trellis,
                                                                   const bool is_closed)
  : Decoder_SIHO<B, Q>(K, is_closed ? 2 * K + 2 * static_cast<int>(std::log2(trellis[0].size())) : 2 * K)
  , m_n_states(static_cast<int>(trellis[0].size()))
  , m_n_memories(static_cast<int>(std::log2(m_n_states)))
  , m_K(K)
  , m_N(is_closed ? 2 * K + 2 * m_n_memories : N)
  , m_L(L)
  , m_n_steps(K + m_n_memories)
  , m_trellis(trellis)
  , m_T({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_T_inv({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_C({ std::vector<int>(m_n_states), std::vector<int>(m_n_states) })
  , m_P(std::vector<Q>(L * m_n_states * (m_n_steps + 1)))
  , m_branch_metr(std::vector<Q>((m_K + m_n_memories) * m_n_states))
  , m_backwards_path(std::vector<int>(L * m_n_states * (m_n_steps + 1)))
  , m_previous_rank(std::vector<int>(L * m_n_states * (m_n_steps + 1)))
  , m_step_result(std::vector<Q>(L * m_n_states))
  , m_M(std::vector<int>(m_n_states * m_n_states))
  , m_decoded(std::vector<std::vector<int>>(L))
  , m_bin_vals({ { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } })
  , m_closing_path(std::vector<int>(m_n_states))
  , m_closing_inputs(std::vector<int>(m_n_states))
  , m_crc(crc.clone())
{
    const std::string name = "Decoder_Viterbi_list_parallel";
    this->set_name(name);
    for (auto& t : this->tasks)
        t->set_replicability(true);

    this->setup();
}

template<typename B, typename Q>
Decoder_Viterbi_list_parallel<B, Q>*
Decoder_Viterbi_list_parallel<B, Q>::clone() const
{
    auto m = new Decoder_Viterbi_list_parallel<B, Q>(*this);
    m->deep_copy(*this);
    m->set_crc(m_crc->clone());
    return m;
}

template<typename B, typename Q>
Decoder_Viterbi_list_parallel<B, Q>::~Decoder_Viterbi_list_parallel()
{
}

/*!
 * \brief Setup all the memory
 */
template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::setup()
{
    std::copy(m_trellis[6].begin(), m_trellis[6].end(), m_T[0].begin());
    std::copy(m_trellis[8].begin(), m_trellis[8].end(), m_T[1].begin());

    std::copy(m_trellis[7].begin(), m_trellis[7].end(), m_C[0].begin());
    std::copy(m_trellis[9].begin(), m_trellis[9].end(), m_C[1].begin());

    int etat_arrivee;
    for (auto i_state_depart = 0; i_state_depart < m_n_states; i_state_depart++)
    {
        for (auto input_bit = 0; input_bit < 2; input_bit++)
        {
            etat_arrivee = m_T[input_bit][i_state_depart];
            m_M[i_state_depart * m_n_states + etat_arrivee] = input_bit + 1;
        }
    }

    for (auto i = 0; i < m_L; i++)
    {
        m_decoded[i] = std::vector<int>(m_K + m_n_memories);
    }

    // Inversion de T
    for (auto bit_sys = 0; bit_sys < 2; bit_sys++)
    {
        for (auto state = 0; state < m_n_states; state++)
        {
            m_T_inv[bit_sys][m_T[bit_sys][state]] = state;
        }
    }

    // Creation des chemins de fermeture
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

    // inversion des chemins de fermeture
    for (auto i_next_state = 0; i_next_state < m_n_states; i_next_state++)
    {
        m_closing_T.push_back(std::vector<int>());
        m_closing_bit_sys.push_back(std::vector<int>());
        for (auto i_prev_state = 0; i_prev_state < m_n_states; i_prev_state++)
        {
            if (m_T[0][i_prev_state] == i_next_state)
            {
                m_closing_T[i_next_state].push_back(i_prev_state);
                m_closing_bit_sys[i_next_state].push_back(0);
            }
            if (m_T[1][i_prev_state] == i_next_state)
            {
                m_closing_T[i_next_state].push_back(i_prev_state);
                m_closing_bit_sys[i_next_state].push_back(1);
            }
        }
    }
}

/*!
 * \brief Reset the matrices
 *
 * \details Set all the values in the node weights' matrice and the backwards
 * path matrice to +inf and the decoded matrice to zeroes.
 */
template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::__reset()
{
    std::fill(m_P.begin(), m_P.end(), gDOUBLE_INF);
    std::fill(m_backwards_path.begin(), m_backwards_path.end(), -1);
    for (auto i = 0; i < m_L; i++)
        std::fill(m_decoded[i].begin(), m_decoded[i].end(), 0);
    std::fill(m_branch_metr.begin(), m_branch_metr.end(), gDOUBLE_INF);
}

// https://stackoverflow.com/a/12399290
// Trie les indices d'un vecteur en fonction de ses valeurs
template<typename T>
std::vector<size_t>
argsort(const std::vector<T>& v)
{

    std::vector<size_t> idx(v.size());
    // Remplit le vecteur de valeurs croissantes, demarrant a 0
    std::iota(idx.begin(), idx.end(), 0);

    // Trie les valeurs de idx.begin() a idx.end() selon la condition
    // specifiee sur les valeurs de v
    std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}

/**
 * @brief Pour le passage en avant, on travaille a etat suivant fixe et on balaie les N etats precedents
 * A chaque etape de decodage, NL nouvelles metriques de noeud sont calculees pour l'etat suivant. Les L
 * meilleures sont conservees et placees dans la liste.
 *
 * @param Y_N LLR du canal
 */
template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::__forward_pass(const Q* Y_N)
{
    const int n_bits = 2; /* juste pour retirer les valeurs hard-codees */
    // initialisation du treillis a l'etat 0 et au rang 0
    m_P[0] = 0.;

    // Decodage "arriere" des K bits d'information du treillis
    for (int i_step = 1; i_step <= m_K; i_step++)
    {
        // Extraction des deux LLR correspondant a l'etape donnee

        std::array<Q, 2> channel_input{ Y_N[2 * (i_step - 1)], Y_N[2 * (i_step - 1) + 1] };

        // reinitialisation des metriques de branche pour cette etape
        std::fill(m_branch_metr.begin(), m_branch_metr.end(), gDOUBLE_INF);

        for (int i_next_state = 0; i_next_state < m_n_states; i_next_state++)
        {

            // Reinitialisation de la matrice step_result
            std::fill(m_step_result.begin(), m_step_result.end(), gDOUBLE_INF);

            for (int bit_sys = 0; bit_sys < n_bits; bit_sys++)
            {
                for (int previous_state = 0; previous_state < m_n_states; previous_state++)
                {
                    for (int i_rank = 0; i_rank < m_L; i_rank++)
                    {
                        const int rank_offset = i_rank * (m_n_steps + 1) * m_n_states;
                        const int idx = rank_offset + (i_step - 1) * m_n_states + previous_state;
                        // Si on peut partir de l'etat donne a l'etape donnee
                        if (m_T[bit_sys][previous_state] == i_next_state && m_P[idx] != gDOUBLE_INF)
                        {
                            // Recuperation du poids du noeud
                            const Q node_weight = m_P[idx];
                            const std::vector<Q> output = m_bin_vals[2 * bit_sys + m_C[bit_sys][previous_state]];
                            Q branch_weight = channel_input[0] * output[0] + channel_input[1] * output[1];
                            ;

                            const Q new_weight = node_weight + branch_weight;

                            m_step_result[i_rank * m_n_states + previous_state] = new_weight;
                        }
                    }
                }
            }
            _process_step(i_step, i_next_state);
        }
    }
}

template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::_process_step(const int i_step, const int i_next_state)
{
    // Tri des valeurs dans step_result
    std::vector<size_t> idx_minima = argsort(m_step_result);

    // Separation des indices en rang et en etats precedents
    for (auto l = 0; l < m_L; l++)
    {
        const int idx =
          l * (m_n_steps + 1) * m_n_states + i_step * m_n_states + i_next_state; // rang ii, etape i_step, etat i_state
        const size_t rank = idx_minima[l] / m_n_states;                          // floor (quotient)
        const size_t prev_state = idx_minima[l] % m_n_states;                    // remainder
        // Si on peut *reduire* la metrique de noeud pour le rang, l'etape et l'etat donne
        if (m_P[idx] > m_step_result[rank * m_n_states + prev_state])
        {
            m_P[idx] = m_step_result[rank * m_n_states + prev_state];
            m_backwards_path[idx] = prev_state;
            m_previous_rank[idx] = rank;
        }
    }
}

template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::__forward_pass_closing(const Q* Y_N)
{
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = FERMETURE
    // Toujours en regardant les etats de depart potentiels depuis l'etat d'arrivee.
    // Cette fois au lieu d'utiliser m_P pour savoir si l'etat precedent etait atteignable,
    // on regarde aussi dans m_closing_T si la transition permet une fermeture de treillis
    // Si c'est le cas, on y fait le traitement adequat
    //
    // Notes a moi meme
    //  Le terme "node" indique un noeud du treillis. Il a pour coordonnees l'etape de decodage
    //  et l'etat correspondant (donc au total, de 0 a (n_steps + 1) * n_states - 1)
    //  Le terme "state" indique un etat du treillis. Il a pour valeurs tout l'intervalle
    //  de 0 a (2 ^ m) - 1 avec m le nombre de memoires.
    //  L'indicateur 'i_' indique un iterateur
    //  L'indicateur 'm_' indique un membre (attribut)
    //  L'indicateur 'g' indique une variable globale
    int n_memories_set_to_0 = 0;

    for (auto i_step = m_K + 1; i_step <= m_K + m_n_memories + 1; i_step++)
    {
        std::array<Q, 2> channel_input{ Y_N[2 * (i_step - 1)], Y_N[2 * (i_step - 1) + 1] };

        std::fill(m_branch_metr.begin(), m_branch_metr.end(), gDOUBLE_INF);
        n_memories_set_to_0++;
        // Pas de saturation ici, on fait confiance au demodulateur pour en appliquer une
        for (auto i_next_state = 0; i_next_state < m_n_states; i_next_state++)
        {
            // Matrice des N*L resultats
            std::fill(m_step_result.begin(), m_step_result.end(), gDOUBLE_INF);
            if (i_next_state <= std::pow(2, m_n_memories - n_memories_set_to_0) - 1)
            {
                const std::vector<int> previous_states = m_closing_T[i_next_state];
                const std::vector<int> inputs = m_closing_bit_sys[i_next_state];
                for (size_t previous_state_idx = 0; previous_state_idx < previous_states.size(); previous_state_idx++)
                {
                    const int previous_state = previous_states[previous_state_idx];
                    const int bit_sys = inputs[previous_state_idx];
                    for (auto i_rank = 0; i_rank < m_L; i_rank++)
                    {
                        const int rank_offset = i_rank * (m_n_steps + 1) * m_n_states;
                        const int prev_node_idx = rank_offset + (i_step - 1) * m_n_states + previous_state;
                        if (m_P[prev_node_idx] != gDOUBLE_INF)
                        {
                            const Q node_weight = m_P[prev_node_idx];

                            Q branch_weight = m_branch_metr[i_next_state * (m_K + m_n_memories) + previous_state];

                            if (branch_weight == gDOUBLE_INF) // pas calcule
                            {
                                const std::vector<Q> output = m_bin_vals[2 * bit_sys + m_C[bit_sys][previous_state]];
                                branch_weight = channel_input[0] * output[0] + channel_input[1] * output[1];
                                m_branch_metr[i_next_state * (m_K + m_n_memories) + previous_state] =
                                  branch_weight; // stockage
                            }

                            const Q new_weight = node_weight + branch_weight;

                            m_step_result[i_rank * m_n_states + previous_state] = new_weight;
                        }
                    }
                }
            }

            _process_step(i_step, i_next_state);
        }
    }
}

template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::_forward_pass(const Q* Y_N)
{
    __forward_pass(Y_N);
    __forward_pass_closing(Y_N);
}

/**
 * @brief Main loop for backwards pass. Called for each rank
 *
 */
template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::__backwards_pass(
  const int l,
  std::vector<B>& data) /* write in data, name and type not defined yet */
{
    int state = 0; // initial state at 0
    int previous_rank = l;

    // std::string output {std::to_string(state)};

    // go through all steps backwards
    for (auto i_step = m_K + m_n_memories; i_step > 0; i_step--)
    {
        const int rank_offset = previous_rank * m_n_states * (m_n_steps + 1);

        const int previous_state = m_backwards_path[rank_offset + i_step * m_n_states + state];
        previous_rank = m_previous_rank[rank_offset + i_step * m_n_states + state];
        // output = std::to_string(previous_state) + " <- " + output;
        //  start storing once the trellis isn't closing
        if (i_step - 1 < m_K) data[i_step - 1] = m_M[previous_state * m_n_states + state] - 1;

        state = previous_state;
    }

    // std::cout << output << std::endl;
}

/**
 * @brief Perform up to L backwards pass depending on if the crc is checked or not
 *
 * This version of the algorithm is for RSC codes. Thus, one can assume it will always
 * start to backtrack from state 0
 *
 * @param Y_N
 * @param X_N
 */
template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::_backwards_pass(B* X_N)
{
    std::vector<B> decoded_path(m_K);
    // compute first path separately, so that we may still return something in the end
    __backwards_pass(0, decoded_path);
    bool stop = m_crc->check(decoded_path);

    std::copy_n(decoded_path.begin(), decoded_path.size(), X_N);
    int l = 1;
    while (!stop && l < m_L)
    {
        __backwards_pass(l, decoded_path);

        stop = m_crc->check(decoded_path);
        // stop = false;
        if (stop) // CRC FOUND!!!
        {
            std::copy_n(decoded_path.begin(), decoded_path.size(), X_N);
        }
        l++;
    }
}

template<typename B, typename Q>
void
Decoder_Viterbi_list_parallel<B, Q>::_backwards_pass(std::vector<B>& X_N)
{
    std::vector<B> decoded_path(m_K);
    // compute first path separately, so that we may still return something in the end
    __backwards_pass(0, decoded_path);
    bool stop = m_crc->check(decoded_path);

    std::copy(decoded_path.begin(), decoded_path.begin() + m_K, X_N.begin());
    int l = 1;
    while (!stop && l < m_L)
    {
        __backwards_pass(l, decoded_path);

        stop = m_crc->check(decoded_path);
        // stop = false;
        if (stop) // CRC FOUND!!!
        {
            std::copy(decoded_path.begin(), decoded_path.begin() + m_K, X_N.begin());
        }
        l++;
    }
}

template<typename B, typename Q>
size_t
Decoder_Viterbi_list_parallel<B, Q>::_argmin(Q* ptr, size_t size)
{
    double minValue = gDOUBLE_INF;
    size_t minIdx = 0;

    for (size_t i = 0; i < size; i++)
    {
        if (minValue > ptr[i])
        {
            minIdx = i;
            minValue = ptr[i];
        }
    }

    return minIdx;
}

template<typename B, typename Q>
size_t
Decoder_Viterbi_list_parallel<B, Q>::_argmin(std::vector<Q> data)
{
    double minValue = gDOUBLE_INF;
    size_t minIdx = 0;

    for (size_t i = 0; i < data.size(); i++)
    {
        if (minValue > data[i])
        {
            minIdx = i;
            minValue = data[i];
        }
    }

    return minIdx;
}

template<typename B, typename Q>
int
Decoder_Viterbi_list_parallel<B, Q>::_decode_siho(const Q* Y_N, B* X_N, const size_t /*frame_id*/)
{
    __reset();
    _forward_pass(Y_N);
    _backwards_pass(X_N);
    return 0;
}

#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Decoder_Viterbi_list_parallel<B_8, Q_8>;
template class aff3ct::module::Decoder_Viterbi_list_parallel<B_16, Q_16>;
template class aff3ct::module::Decoder_Viterbi_list_parallel<B_32, Q_32>;
template class aff3ct::module::Decoder_Viterbi_list_parallel<B_64, Q_64>;
#else
template class aff3ct::module::Decoder_Viterbi_list_parallel<B, Q>;
#endif