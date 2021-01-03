#include <chrono>
#include <iostream>
#include <random>
#include <vector>

// needs -fno-math-errno, so std::sqrt() can be vectorized

#if defined __GNUC__
#    define IVDEP _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
#    define IVDEP __pragma(loop(ivdep))
#elif defined __clang__
#    define IVDEP _Pragma("clang loop vectorize(enable) interleave(enable) distribute(enable)")
#else
#    error "Please define IVDEP for your compiler!"
#endif

struct Stopwatch
{
    using clock = std::chrono::high_resolution_clock;

    auto elapsedAndReset() -> double
    {
        const auto now = clock::now();
        const auto seconds = std::chrono::duration<double>{now - last}.count();
        last = now;
        return seconds;
    }

private:
    clock::time_point last = clock::now();
};

constexpr auto PROBLEM_SIZE = 16 * 1024;
constexpr auto STEPS = 5;
constexpr auto TIMESTEP = 0.0001f;
constexpr auto EPS2 = 0.01f;

namespace AoS
{
    struct Vec
    {
        float x;
        float y;
        float z;

        auto operator*=(float s) -> Vec&
        {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }

        auto operator*=(Vec v) -> Vec&
        {
            x *= v.x;
            y *= v.y;
            z *= v.z;
            return *this;
        }

        auto operator+=(Vec v) -> Vec&
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        auto operator-=(Vec v) -> Vec&
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        friend auto operator+(Vec a, Vec b) -> Vec
        {
            return a += b;
        }

        friend auto operator-(Vec a, Vec b) -> Vec
        {
            return a -= b;
        }

        friend auto operator*(Vec a, float s) -> Vec
        {
            return a *= s;
        }

        friend auto operator*(Vec a, Vec b) -> Vec
        {
            return a *= b;
        }
    };

    struct Particle
    {
        Vec pos;
        Vec vel;
        float mass;
    };

    inline void pPInteraction(Particle& pi, const Particle& pj)
    {
        const Vec distance = pi.pos - pj.pos;
        const Vec distanceSqr = distance * distance;
        const float distSqr = EPS2 + distanceSqr.x + distanceSqr.y + distanceSqr.z;
        const float distSixth = distSqr * distSqr * distSqr;
        const float invDistCube = 1.0f / std::sqrt(distSixth);
        const float sts = pj.mass * invDistCube * TIMESTEP;
        pi.vel += distanceSqr * sts;
    }

    void update(Particle* particles)
    {
        IVDEP
        for (std::size_t i = 0; i < PROBLEM_SIZE; i++)
            for (std::size_t j = 0; j < PROBLEM_SIZE; j++)
                pPInteraction(particles[i], particles[j]);
    }

    void move(Particle* particles)
    {
        IVDEP
        for (std::size_t i = 0; i < PROBLEM_SIZE; i++)
            particles[i].pos += particles[i].vel * TIMESTEP;
    }

    void run()
    {
        std::vector<Particle> particles(PROBLEM_SIZE);

        std::default_random_engine engine;
        std::normal_distribution<float> dist(0.0f, 1.0f);
        for (auto& p : particles)
        {
            p.pos.x = dist(engine);
            p.pos.y = dist(engine);
            p.pos.z = dist(engine);
            p.vel.x = dist(engine) / 10.0f;
            p.vel.y = dist(engine) / 10.0f;
            p.vel.z = dist(engine) / 10.0f;
            p.mass = dist(engine) / 100.0f;
        }

        Stopwatch watch;
        double sumUpdate = 0;
        double sumMove = 0;
        for (std::size_t s = 0; s < STEPS; ++s)
        {
            update(particles.data());
            sumUpdate += watch.elapsedAndReset();
            move(particles.data());
            sumMove += watch.elapsedAndReset();
        }
        std::cout << "AoS\t" << sumUpdate / STEPS << '\t' << sumMove / STEPS << '\n';
    }
} // namespace AoS

namespace SoA
{
    template <typename T, std::size_t Alignment>
    struct AlignedAllocator
    {
        using value_type = T;

        inline AlignedAllocator() noexcept = default;

        template <typename T2>
        inline AlignedAllocator(AlignedAllocator<T2, Alignment> const&) noexcept
        {
        }

        inline ~AlignedAllocator() noexcept = default;

        inline auto allocate(std::size_t n) -> T*
        {
            return static_cast<T*>(::operator new[](n * sizeof(T), std::align_val_t{Alignment}));
        }

        inline void deallocate(T* p, std::size_t)
        {
            ::operator delete[](p, std::align_val_t{Alignment});
        }

        template <typename T2>
        struct rebind
        {
            using other = AlignedAllocator<T2, Alignment>;
        };

        auto operator!=(const AlignedAllocator<T, Alignment>& other) const -> bool
        {
            return !(*this == other);
        }

        auto operator==(const AlignedAllocator<T, Alignment>& other) const -> bool
        {
            return true;
        }
    };

    inline void pPInteraction(
        float piposx,
        float piposy,
        float piposz,
        float& pivelx,
        float& pively,
        float& pivelz,
        float pjposx,
        float pjposy,
        float pjposz,
        float pjmass)
    {
        const float xdistance = piposx - pjposx;
        const float ydistance = piposy - pjposy;
        const float zdistance = piposz - pjposz;
        const float xdistanceSqr = xdistance * xdistance;
        const float ydistanceSqr = ydistance * ydistance;
        const float zdistanceSqr = zdistance * zdistance;
        const float distSqr = EPS2 + xdistanceSqr + ydistanceSqr + zdistanceSqr;
        const float distSixth = distSqr * distSqr * distSqr;
        const float invDistCube = 1.0f / std::sqrt(distSixth);
        const float sts = pjmass * invDistCube * TIMESTEP;
        pivelx += xdistanceSqr * sts;
        pively += ydistanceSqr * sts;
        pivelz += zdistanceSqr * sts;
    }

    void update(
        const float* posx,
        const float* posy,
        const float* posz,
        float* velx,
        float* vely,
        float* velz,
        const float* mass)
    {
        IVDEP
        for (std::size_t i = 0; i < PROBLEM_SIZE; i++)
            for (std::size_t j = 0; j < PROBLEM_SIZE; j++)
                pPInteraction(posx[i], posy[i], posz[i], velx[i], vely[i], velz[i], posx[j], posy[j], posz[j], mass[j]);
    }

    void move(float* posx, float* posy, float* posz, const float* velx, const float* vely, const float* velz)
    {
        IVDEP
        for (std::size_t i = 0; i < PROBLEM_SIZE; i++)
        {
            posx[i] += velx[i] * TIMESTEP;
            posy[i] += vely[i] * TIMESTEP;
            posz[i] += velz[i] * TIMESTEP;
        }
    }

    void run()
    {
        using Vector = std::vector<float, AlignedAllocator<float, 64>>;
        Vector posx(PROBLEM_SIZE);
        Vector posy(PROBLEM_SIZE);
        Vector posz(PROBLEM_SIZE);
        Vector velx(PROBLEM_SIZE);
        Vector vely(PROBLEM_SIZE);
        Vector velz(PROBLEM_SIZE);
        Vector mass(PROBLEM_SIZE);

        std::default_random_engine engine;
        std::normal_distribution<float> dist(0.0f, 1.0f);
        for (std::size_t i = 0; i < PROBLEM_SIZE; ++i)
        {
            posx[i] = dist(engine);
            posy[i] = dist(engine);
            posz[i] = dist(engine);
            velx[i] = dist(engine) / 10.0f;
            vely[i] = dist(engine) / 10.0f;
            velz[i] = dist(engine) / 10.0f;
            mass[i] = dist(engine) / 100.0f;
        }

        Stopwatch watch;
        double sumUpdate = 0;
        double sumMove = 0;
        for (std::size_t s = 0; s < STEPS; ++s)
        {
            update(posx.data(), posy.data(), posz.data(), velx.data(), vely.data(), velz.data(), mass.data());
            sumUpdate += watch.elapsedAndReset();
            move(posx.data(), posy.data(), posz.data(), velx.data(), vely.data(), velz.data());
            sumMove += watch.elapsedAndReset();
        }
        std::cout << "SoA\t" << sumUpdate / STEPS << '\t' << sumMove / STEPS << '\n';
    }
} // namespace SoA

namespace AoSoA
{
    constexpr auto LANES = 16;
    constexpr auto BLOCKS = PROBLEM_SIZE / LANES;

    struct alignas(64) ParticleBlock
    {
        struct
        {
            float x[LANES];
            float y[LANES];
            float z[LANES];
        } pos;
        struct
        {
            float x[LANES];
            float y[LANES];
            float z[LANES];
        } vel;
        float mass[LANES];
    };

    inline void pPInteraction(
        float piposx,
        float piposy,
        float piposz,
        float& pivelx,
        float& pively,
        float& pivelz,
        float pjposx,
        float pjposy,
        float pjposz,
        float pjmass)
    {
        auto xdistance = piposx - pjposx;
        auto ydistance = piposy - pjposy;
        auto zdistance = piposz - pjposz;
        xdistance *= xdistance;
        ydistance *= ydistance;
        zdistance *= zdistance;
        const float distSqr = EPS2 + xdistance + ydistance + zdistance;
        const float distSixth = distSqr * distSqr * distSqr;
        const float invDistCube = 1.0f / std::sqrt(distSixth);
        const float sts = pjmass * invDistCube * TIMESTEP;
        pivelx += xdistance * sts;
        pively += ydistance * sts;
        pivelz += zdistance * sts;
    }

    void update(ParticleBlock* particles)
    {
        for (std::size_t bi = 0; bi < BLOCKS; bi++)
            for (std::size_t bj = 0; bj < BLOCKS; bj++)
                for (std::size_t j = 0; j < LANES; j++)
                {
                    IVDEP
                    for (std::size_t i = 0; i < LANES; i++)
                    {
                        auto& blockI = particles[bi];
                        const auto& blockJ = particles[bj];
                        pPInteraction(
                            blockI.pos.x[i],
                            blockI.pos.y[i],
                            blockI.pos.z[i],
                            blockI.vel.x[i],
                            blockI.vel.y[i],
                            blockI.vel.z[i],
                            blockJ.pos.x[j],
                            blockJ.pos.y[j],
                            blockJ.pos.z[j],
                            blockJ.mass[j]);
                    }
                }
    }

    void move(ParticleBlock* particles)
    {
        for (std::size_t bi = 0; bi < BLOCKS; bi++)
        {
            IVDEP
            for (std::size_t i = 0; i < LANES; ++i)
            {
                auto& block = particles[bi];
                block.pos.x[i] += block.vel.x[i] * TIMESTEP;
                block.pos.y[i] += block.vel.y[i] * TIMESTEP;
                block.pos.z[i] += block.vel.z[i] * TIMESTEP;
            }
        }
    }

    void run()
    {
        std::vector<ParticleBlock> particles(BLOCKS);

        std::default_random_engine engine;
        std::normal_distribution<float> dist(0.0f, 1.0f);
        for (std::size_t bi = 0; bi < BLOCKS; ++bi)
        {
            auto& block = particles[bi];
            for (std::size_t i = 0; i < LANES; ++i)
            {
                block.pos.x[i] = dist(engine);
                block.pos.y[i] = dist(engine);
                block.pos.z[i] = dist(engine);
                block.vel.x[i] = dist(engine) / 10.0f;
                block.vel.y[i] = dist(engine) / 10.0f;
                block.vel.z[i] = dist(engine) / 10.0f;
                block.mass[i] = dist(engine) / 100.0f;
            }
        }

        Stopwatch watch;
        double sumUpdate = 0;
        double sumMove = 0;
        for (std::size_t s = 0; s < STEPS; ++s)
        {
            update(particles.data());
            sumUpdate += watch.elapsedAndReset();
            move(particles.data());
            sumMove += watch.elapsedAndReset();
        }
        std::cout << "AoSoA\t" << sumUpdate / STEPS << '\t' << sumMove / STEPS << '\n';
    }
} // namespace AoSoA

int main()
{
    std::cout << PROBLEM_SIZE / 1000 << "k particles "
              << "(" << PROBLEM_SIZE * sizeof(float) * 7 / 1024 << "kiB)\n";

    AoS::run();
    SoA::run();
    AoSoA::run();

    return 0;
}
