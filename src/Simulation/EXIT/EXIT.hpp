#ifndef SIMULATION_EXIT_HPP_
#define SIMULATION_EXIT_HPP_

#include <vector>
#include <mipp.h>

#include "Module/Source/Source.hpp"
#include "Module/Encoder/Encoder.hpp"
#include "Module/Modem/Modem.hpp"
#include "Module/Channel/Channel.hpp"
#include "Module/Decoder/Decoder_SISO.hpp"

#include "Tools/Display/Terminal/EXIT/Terminal_EXIT.hpp"
#include "Tools/Codec/Codec_SISO.hpp"

#include "Factory/Simulation/EXIT/EXIT.hpp"

#include "../Simulation.hpp"

namespace aff3ct
{
namespace simulation
{
template <typename B = int, typename R = float>
class EXIT : public Simulation
{
protected:
	const factory::EXIT::parameters &params; // simulation parameters
	tools::Codec_SISO<B,R> &codec;

	// channel gains
	mipp::vector<R> H_N;

	// data vectors
	mipp::vector<B> B_K, B_N, X_N1;
	mipp::vector<R> X_K, X_N2;
	mipp::vector<R> Y_N, Y_K;
	mipp::vector<R> La_K1;
	mipp::vector<R> Lch_N1;
	mipp::vector<R> La_K2;
	mipp::vector<R> Lch_N2;
	mipp::vector<R> Le_K;
	mipp::vector<R> sys, par;
	mipp::vector<B> B_buff;
	mipp::vector<R> Le_buff, La_buff;

	// EXIT simu parameters
	const int n_trials;
	int cur_trial;
	double I_A, I_E;

	// code specifications
	float sig_a;
	float sigma;
	float snr;

	// communication chain
	module::Source<B>       *source;
	module::Encoder<B>      *encoder;
	module::Modem<B,R,R>    *modem;
	module::Modem<B,R,R>    *modem_a;
	module::Channel<R>      *channel;
	module::Channel<R>      *channel_a;
	module::Decoder_SISO<R> *siso;
	tools::Terminal_EXIT    *terminal;

public:
	EXIT(const factory::EXIT::parameters& params, tools::Codec_SISO<B,R> &codec);
	virtual ~EXIT();
	
	void launch();

private:
	void build_communication_chain();
	void simulation_loop          ();

	static double measure_mutual_info_avg  (const mipp::vector<R>& llrs, const mipp::vector<B>& bits);
	static double measure_mutual_info_histo(const mipp::vector<R>& llrs, const mipp::vector<B>& bits);

protected:
	virtual void release_objects  ();

	virtual module::Source<B>*       build_source   (              );
	virtual module::Encoder<B>*      build_encoder  (              );
	virtual module::Modem  <B,R,R>*  build_modem    (              );
	virtual module::Modem  <B,R,R>*  build_modem_a  (              );
	virtual module::Channel<R>*      build_channel  (const int size);
	virtual module::Channel<R>*      build_channel_a(const int size);
	virtual module::Decoder_SISO<R>* build_siso     (              );
	        tools::Terminal_EXIT*    build_terminal (              );
};
}
}

#endif /* SIMULATION_EXIT_HPP_ */