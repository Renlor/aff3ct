/*!
 * \file
 * \brief Class factory::Codec_RSC.
 */
#ifndef FACTORY_CODEC_RSC_HPP
#define FACTORY_CODEC_RSC_HPP

#include <string>
#include <map>
#include <cli.hpp>

#include "Module/CRC/CRC.hpp"
#include "Module/Codec/RSC/Codec_RSC.hpp"
#include "Factory/Module/Codec/Codec_SISO_SIHO.hpp"

namespace aff3ct
{
namespace factory
{
extern const std::string Codec_RSC_name;
extern const std::string Codec_RSC_prefix;
class Codec_RSC : public Codec_SISO_SIHO
{
public:
	explicit Codec_RSC(const std::string &p = Codec_RSC_prefix);
	virtual ~Codec_RSC() = default;
	Codec_RSC* clone() const;

	// parameters construction
	void get_description(cli::Argument_map_info &args) const;
	void store          (const cli::Argument_map_value &vals);
	void get_headers    (std::map<std::string,header_list>& headers, const bool full = true) const;

	// builder
	template <typename B = int, typename Q = float>
	module::Codec_RSC<B,Q>* build(module::CRC<B>* crc = nullptr) const;
};
}
}

#endif /* FACTORY_CODEC_RSC_HPP */
