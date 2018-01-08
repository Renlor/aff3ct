#ifndef ARGUMENT_TYPE_INTEGER_HPP_
#define ARGUMENT_TYPE_INTEGER_HPP_

#include <string>
#include <stdexcept>

#include "Argument_type.hpp"

namespace aff3ct
{
namespace tools
{

template <typename T = int, typename... Ranges>
class Integer_type : public Argument_type_limited<T,Ranges...>
{
public:
	Integer_type(const Ranges*... ranges)
	: Argument_type_limited<T,Ranges...>("integer", ranges...)
	{ }

	virtual ~Integer_type() {};

	virtual Integer_type<T,Ranges...>* clone() const
	{
		auto clone = new Integer_type<T,Ranges...>(*this);

		return dynamic_cast<Integer_type<T,Ranges...>*>(this->clone_ranges(clone));
	}

	virtual T convert(const std::string& val) const
	{
		return (T)std::stoi(val);
	}

	virtual void check(const std::string& val) const
	{
		T int_num;

		try
		{
			int_num = this->convert(val);
		}
		catch(std::exception&)
		{
			throw std::runtime_error("shall be an integer");
		}

		this->check_ranges(int_num);
	}
};

template <typename T = int, typename... Ranges>
Integer_type<T,Ranges...>* Integer(Ranges*... ranges)
{
	return new Integer_type<T,Ranges...>(ranges...);
}

}
}
#endif /* ARGUMENT_TYPE_INTEGER_HPP_ */