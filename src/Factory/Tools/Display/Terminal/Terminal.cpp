#include <streampu.hpp>
#include <utility>

#include "Factory/Tools/Display/Terminal/Terminal.hpp"
#include "Tools/Display/Terminal/Standard/Terminal_std.hpp"
#include "Tools/Documentation/documentation.h"

using namespace aff3ct;
using namespace aff3ct::factory;

const std::string aff3ct::factory::Terminal_name = "Terminal";
const std::string aff3ct::factory::Terminal_prefix = "ter";

Terminal ::Terminal(const std::string& prefix)
  : Factory(Terminal_name, Terminal_name, prefix)
{
}

Terminal ::Terminal(const std::string& name, const std::string& prefix)
  : Factory(name, Terminal_name, prefix)
{
}

Terminal*
Terminal ::clone() const
{
    return new Terminal(*this);
}

void
Terminal ::get_description(cli::Argument_map_info& args) const
{
    auto p = this->get_prefix();
    const std::string class_name = "factory::Terminal::";

    tools::add_arg(args, p, class_name + "p+type", cli::Text(cli::Including_set("STD")));

    tools::add_arg(args, p, class_name + "p+no", cli::None());

    tools::add_arg(args, p, class_name + "p+freq", cli::Integer(cli::Positive()));
}

void
Terminal ::store(const cli::Argument_map_value& vals)
{
    auto p = this->get_prefix();

    if (vals.exist({ p + "-type" })) this->type = vals.at({ p + "-type" });
    if (vals.exist({ p + "-no" })) this->disabled = true;
    if (vals.exist({ p + "-freq" })) this->frequency = std::chrono::milliseconds(vals.to_int({ p + "-freq" }));
}

void
Terminal ::get_headers(std::map<std::string, tools::header_list>& headers, const bool /*full*/) const
{
    auto p = this->get_prefix();

    headers[p].push_back(std::make_pair("Enabled", this->disabled ? "no" : "yes"));
    headers[p].push_back(std::make_pair("Frequency (ms)", std::to_string(this->frequency.count())));
}

spu::tools::Terminal*
Terminal ::build(const std::vector<std::unique_ptr<spu::tools::Reporter>>& reporters) const
{
    if (this->type == "STD") return new spu::tools::Terminal_std(reporters);

    throw spu::tools::cannot_allocate(__FILE__, __LINE__, __func__);
}