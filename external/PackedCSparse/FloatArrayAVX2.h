template<size_t k>
        struct m256dArray
{
    __m256d x[k];
    
    m256dArray() {};
    
    m256dArray(const double rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = _mm256_set1_pd(rhs);
    }
    
    template<size_t k2>
            m256dArray(const m256dArray<k2>& rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = rhs.x[i % k2];
    }
    
    m256dArray operator+(const m256dArray& rhs) const
    {
        m256dArray out;
        for (size_t i = 0; i < k; i++)
            out.x[i] = _mm256_add_pd(x[i], rhs.x[i]);
        return out;
    }
    
    m256dArray operator-(const m256dArray& rhs) const
    {
        m256dArray out;
        for (size_t i = 0; i < k; i++)
            out.x[i] = _mm256_sub_pd(x[i], rhs.x[i]);
        return out;
    }
    
    m256dArray operator*(const m256dArray& rhs) const
    {
        m256dArray out;
        for (size_t i = 0; i < k; i++)
            out.x[i] = _mm256_mul_pd(x[i], rhs.x[i]);
        return out;
    }
    
    m256dArray operator/(const m256dArray& rhs) const
    {
        m256dArray out;
        for (size_t i = 0; i < k; i++)
            out.x[i] = _mm256_div_pd(x[i], rhs.x[i]);
        return out;
    }
    
    m256dArray& operator+=(const m256dArray& rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = _mm256_add_pd(x[i], rhs.x[i]);
        return *this;
    }
    
    m256dArray& operator-=(const m256dArray& rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = _mm256_sub_pd(x[i], rhs.x[i]);
        return *this;
    }
    
    m256dArray& operator*=(const m256dArray& rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = _mm256_mul_pd(x[i], rhs.x[i]);
        return *this;
    }
    
    m256dArray& operator/=(const m256dArray& rhs)
    {
        for (size_t i = 0; i < k; i++)
            x[i] = _mm256_div_pd(x[i], rhs.x[i]);
        return *this;
    }
    
    explicit operator bool() const
    {
        bool ret = false;
        __m256d z = _mm256_set1_pd(0.0);
        for (size_t i = 0; i < k; i++)
        {
            __m256d c = _mm256_cmp_pd(x[i], z, _CMP_EQ_OQ);
            ret = ret || (_mm256_movemask_pd(c) != 0xf);
        }
        return ret;
    }
    
    static double get(const m256dArray& x, size_t index)
    {
        double y[4];
        _mm256_store_pd(y, x.x[index / 4]);
        return y[index & 3];
    }
    
    static void set(m256dArray& x, size_t index, double value)
    {
        __m256d v = _mm256_broadcast_sd(&value);
        switch (index & 3)
        {
            case 0:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 1); break;
            case 1:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 2); break;
            case 2:  x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 4); break;
            default: x.x[index / 4] = _mm256_blend_pd(x.x[index / 4], v, 8); break;
        }
    }
    
    static m256dArray abs(const m256dArray& x)
    {
        const __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
        
        m256dArray out;
        for (size_t i = 0; i < k; i++)
            out.x[i] = _mm256_and_pd(x.x[i], mask);
        return out;
    }
    
    static m256dArray log(const m256dArray& x)
    {
        // gcc does not support _mm256_log_pd
        // Do it sequentially instead
        
        //m256dArray out;
        //for (size_t i = 0; i < k; i++)
        //	out.x[i] = _mm256_log_pd(x.x[i]);
        
        m256dArray out;
        for (size_t i = 0; i < 4*k; i++)
            set(out, i, std::log(get(x,i)));
        return out;
    }
    
    static void fmadd(m256dArray& a, const m256dArray& b, const double& c)
    {
        auto cx = _mm256_set1_pd(c);
        for (size_t i = 0; i < k; i++)
            a.x[i] = _mm256_fmadd_pd(b.x[i], cx, a.x[i]);
    }
    
    static void fnmadd(m256dArray& a, const m256dArray& b, const double& c)
    {
        auto cx = _mm256_set1_pd(c);
        for (size_t i = 0; i < k; i++)
            a.x[i] = _mm256_fnmadd_pd(b.x[i], cx, a.x[i]);
    }
    
    static void fmadd(m256dArray& a, const m256dArray& b, const m256dArray& c)
    {
        for (size_t i = 0; i < k; i++)
            a.x[i] = _mm256_fmadd_pd(b.x[i], c.x[i], a.x[i]);
    }
    
    static void fnmadd(m256dArray& a, const m256dArray& b, const m256dArray& c)
    {
        for (size_t i = 0; i < k; i++)
            a.x[i] = _mm256_fnmadd_pd(b.x[i], c.x[i], a.x[i]);
    }
    
    static m256dArray clipped_sqrt(const m256dArray& x, const double nonpos_output)
    {
        m256dArray out;
        
        const __m256d large = { nonpos_output, nonpos_output, nonpos_output, nonpos_output };
        const __m256d zero = _mm256_setzero_pd();
        for (size_t i = 0; i < k; i++)
        {
            __m256d xi = x.x[i];
            __m256d mask = _mm256_cmp_pd(xi, zero, _CMP_LE_OS); // mask = (rhs.x[i]<= 0) ? -1 : 0
            out.x[i] = _mm256_blendv_pd(_mm256_sqrt_pd(xi), large, mask);
        }
        return out;
    }
    
    static m256dArray sign(std::mt19937_64& gen)
    {
        m256dArray out;
        const __m256i bits = _mm256_set_epi64x(1, 2, 4, 8);
        const __m256d zero = _mm256_setzero_pd();
        const __m256d pos = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
        const __m256d neg = _mm256_set_pd(-1.0, -1.0, -1.0, -1.0);
        
        unsigned long long seed = gen();
        for (size_t i = 0; i < k; i++)
        {
            __m256i s = _mm256_set1_epi64x((seed >> (4 * i)) & 15);
            __m256i xi = _mm256_and_si256(s, bits);
            __m256d x = _mm256_castsi256_pd(xi);
            __m256d mask = _mm256_cmp_pd(x, zero, _CMP_EQ_OQ); // mask = (rhs.x[i] == 0) ? -1 : 0
            out.x[i] = _mm256_blendv_pd(pos, neg, mask);
            if ((i & 63) == 63) seed = gen();
        }
        return out;
    }
};

template <size_t k>
        struct FloatTypeSelector<double, k>
{
    static_assert(k == 1 || k % 4 == 0, "Array<double,k> assumes k = 1 or a multiple of 4");
    using type = typename std::conditional< k == 1, double, m256dArray<k / 4>>::type;
    using funcImpl = typename std::conditional< k == 1, BaseScalarImpl<double>, m256dArray<k / 4>>::type;
};

template <size_t k, size_t l>
        struct FloatTypeSelector<m256dArray<k>, l>
{
    using type = m256dArray<k* l>;
    using funcImpl = m256dArray<k* l>;
};
