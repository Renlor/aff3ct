#include <sstream>
#include <string>

#include "Tools/Exception/exception.hpp"
#include "Tools/Math/utils.h"
#include "Factory/Module/Puncturer/Puncturer.hpp"
#include "Factory/Module/Encoder/Encoder.hpp"
#include "Module/Codec/RS/Codec_RS.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B, typename Q>
Codec_RS<B,Q>
::Codec_RS(const factory::Encoder_RS &enc_params,
           const factory::Decoder_RS &dec_params)
: Codec          <B,Q>(enc_params.K, enc_params.N_cw, enc_params.N_cw, enc_params.tail_length, enc_params.n_frames),
  Codec_SIHO_HIHO<B,Q>(enc_params.K, enc_params.N_cw, enc_params.N_cw, enc_params.tail_length, enc_params.n_frames),
  GF_poly(tools::next_power_of_2(dec_params.N_cw) -1, dec_params.t)
{
	const std::string name = "Codec_RS";
	this->set_name(name);

	// ----------------------------------------------------------------------------------------------------- exceptions
	if (enc_params.K != dec_params.K)
	{
		std::stringstream message;
		message << "'enc_params.K' has to be equal to 'dec_params.K' ('enc_params.K' = " << enc_params.K
		        << ", 'dec_params.K' = " << dec_params.K << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (enc_params.N_cw != dec_params.N_cw)
	{
		std::stringstream message;
		message << "'enc_params.N_cw' has to be equal to 'dec_params.N_cw' ('enc_params.N_cw' = " << enc_params.N_cw
		        << ", 'dec_params.N_cw' = " << dec_params.N_cw << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (enc_params.n_frames != dec_params.n_frames)
	{
		std::stringstream message;
		message << "'enc_params.n_frames' has to be equal to 'dec_params.n_frames' ('enc_params.n_frames' = "
		        << enc_params.n_frames << ", 'dec_params.n_frames' = " << dec_params.n_frames << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	// ---------------------------------------------------------------------------------------------------- allocations
	factory::Puncturer pct_params;
	pct_params.type     = "NO";
	pct_params.K        = enc_params.K    * dec_params.m;
	pct_params.N        = enc_params.N_cw * dec_params.m;
	pct_params.N_cw     = enc_params.N_cw * dec_params.m;
	pct_params.n_frames = enc_params.n_frames;

	this->set_puncturer(pct_params.build<B,Q>());
	try
	{
		this->set_encoder(enc_params.build<B>(GF_poly));
	}
	catch (tools::cannot_allocate const&)
	{
		this->set_encoder(static_cast<const factory::Encoder*>(&enc_params)->build<B>());
	}
	if (dec_params.implem == "GENIUS")
		this->get_encoder()->set_memorizing(true);
	this->set_decoder_siho_hiho(dec_params.build_hiho<B,Q>(GF_poly, this->get_encoder()));
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Codec_RS<B_8,Q_8>;
template class aff3ct::module::Codec_RS<B_16,Q_16>;
template class aff3ct::module::Codec_RS<B_32,Q_32>;
template class aff3ct::module::Codec_RS<B_64,Q_64>;
#else
template class aff3ct::module::Codec_RS<B,Q>;
#endif
// ==================================================================================== explicit template instantiation
