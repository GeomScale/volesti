// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis
// Copyright (c) 2022 Ioannis Iakovidis

// This file is converted from PolytopeSamplerMatlab
//(https://github.com/ConstrainedSampler/PolytopeSamplerMatlab/blob/master/code/solver/PackedCSparse/PackedChol.h) by Ioannis Iakovidis

#pragma once
#include <immintrin.h>
#include <random>
#include <type_traits>
namespace PackedCSparse {
	template <typename T, size_t k>
	struct BaseImpl
	{
		T x[k];

		BaseImpl() {};

		BaseImpl(const T& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] = rhs;
		}

		BaseImpl operator+(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] + rhs.x[i];
			return lhs;
		}

		BaseImpl operator-(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] - rhs.x[i];
			return lhs;
		}

		BaseImpl operator*(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] * rhs.x[i];
			return lhs;
		}

		BaseImpl operator/(const BaseImpl& rhs) const
		{
			BaseImpl lhs;
			for (size_t i = 0; i < k; i++)
				lhs.x[i] = x[i] / rhs.x[i];
			return lhs;
		}

		BaseImpl& operator+=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] += rhs.x[i];
			return *this;
		}

		BaseImpl& operator-=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] -= rhs.x[i];
			return *this;
		}

		BaseImpl& operator*=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] *= rhs.x[i];
			return *this;
		}

		BaseImpl& operator/=(const BaseImpl& rhs)
		{
			for (size_t i = 0; i < k; i++)
				x[i] /= rhs.x[i];
			return *this;
		}

		explicit operator bool() const
		{
			bool ret = false;
			for (size_t i = 0; i < k; i++)
				ret = ret || bool(x[i]);
			return ret;
		}

		static T get(const BaseImpl& a, size_t index)
		{
			return a.x[index];
		}

		static void set(BaseImpl& a, size_t index, const T& value)
		{
			a.x[index] = value;
		}

		static BaseImpl abs(const BaseImpl& a)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
                out.x[i] = std::abs(a.x[i]);
			return out;
		}

		static BaseImpl log(const BaseImpl& a)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
				out.x[i] = std::log(a.x[i]);
			return out;
		}

		static void fmadd(BaseImpl& a, const BaseImpl& b, const BaseImpl& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] += b.x[i] * c.x[i];
		}

		static void fnmadd(BaseImpl& a, const BaseImpl& b, const BaseImpl& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] -= b.x[i] * c.x[i];
		}

		static void fmadd(BaseImpl& a, const BaseImpl& b, const T& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] += b.x[i] * c;
		}

		static void fnmadd(BaseImpl& a, const BaseImpl& b, const T& c)
		{
			for (size_t i = 0; i < k; i++)
				a.x[i] -= b.x[i] * c;
		}

		static BaseImpl clipped_sqrt(const BaseImpl& a, const T nonpos_output)
		{
			BaseImpl out;
			for (size_t i = 0; i < k; i++)
			{
				T r = a.x[i];
				if (r > 0)
					out.x[i] = sqrt(r);
				else
					out.x[i] = nonpos_output;
			}
			return out;
		}

		static BaseImpl sign(std::mt19937_64& gen)
		{
			BaseImpl out;
			unsigned long long seed = gen();
			for (size_t i = 0; i < k; i++)
			{
				out.x[i] = T((2 * ((seed >> i) & 1)) - 1.0);
				if ((i & 63) == 63) seed = gen();
			}
			return out;
		}

	};

	template <typename T>
	struct BaseScalarImpl
	{
		static T get(const T& x, size_t index)
		{
			return x;
		}

		static void set(T& x, size_t index, T& value)
		{
			x = value;
		}

		static T abs(const T &x)
		{
			return ::abs(x);
		}

		static T log(const T &x)
		{
			return ::log(x);
		}

		static void fmadd(T& a, const T& b, const T& c)
		{
			a += b * c;
		}

		static void fnmadd(T& a, const T& b, const T& c)
		{
			a -= b * c;
		}

		static T clipped_sqrt(const T& x, const T& nonpos_output)
		{
			if (x > 0.0)
				return sqrt(x);
			else
				return nonpos_output;
		}

		static T sign(std::mt19937_64& gen)
		{
			unsigned long long seed = gen();
			return T((2 * (seed & 1)) - 1.0);
		}
	};

	template <typename T, size_t k>
	struct FloatTypeSelector
	{
		using type = typename std::conditional<k == 1, T, BaseImpl<T, k>>::type;
		using funcImpl = typename std::conditional<k == 1, BaseScalarImpl<T>, BaseImpl<T, k>>::type;
	};

    #ifdef __AVX2__
        #include "FloatArrayAVX2.h"
    #else
    template <size_t k>
            struct FloatTypeSelector<double, k>
    {
        using type = typename std::conditional< k == 1, double, BaseImpl<double, k>>::type;
        using funcImpl = typename std::conditional< k == 1, BaseScalarImpl<double>, BaseImpl<double, k>>::type;
    };

    template <size_t k, size_t l>
        struct FloatTypeSelector<BaseImpl<double, k>, l>
    {
        using type = BaseImpl<double, k*l>;
        using funcImpl = BaseImpl<double, k*l>;
    };
    #endif

	template <typename T, size_t k, size_t l>
	struct FloatTypeSelector<BaseImpl<T, k>, l>
	{
		using type = BaseImpl<T, k* l>;
		using funcImpl = BaseImpl<T, k* l>;
	};

	template<typename T, size_t k = 1>
	using FloatArray = typename FloatTypeSelector<T, k>::type;

	template<typename T, size_t k = 1>
	using FloatArrayFunc = typename FloatTypeSelector<T, k>::funcImpl;

	template<typename T>
	auto get(const T& a, size_t index) -> decltype(FloatArrayFunc<T>::get(a, index))
	{
		return FloatArrayFunc<T>::get(a, index);
	}

	template<typename T1, typename T2>
	void set(T1& a, size_t index, T2 value)
	{
		FloatArrayFunc<T1>::set(a, index, value);
	}

	template<typename T1, typename T2, typename T3>
	void fmadd(T1& a, const T2& b, const T3& c)
	{
		FloatArrayFunc<T1>::fmadd(a, b, c);
	}

	template<typename T1, typename T2, typename T3>
	void fnmadd(T1& a, const T2& b, const T3& c)
	{
		FloatArrayFunc<T1>::fnmadd(a, b, c);
	}

	template<typename T1, typename T2>
	T1 clipped_sqrt(const T1& a, const T2 b)
	{
		return FloatArrayFunc<T1>::clipped_sqrt(a, b);
	}

	template<typename T>
	T abs(const T& a)
	{
		return FloatArrayFunc<T>::abs(a);
	}

	template<typename T>
	T log(const T& a)
	{
		return FloatArrayFunc<T>::log(a);
	}

	template<typename T>
	T sign(std::mt19937_64& gen)
	{
		return FloatArrayFunc<T>::sign(gen);
	}
}
