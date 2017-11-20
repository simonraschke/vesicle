#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)


#include <vector>
#include <memory>
#define PARTICLERANGE std::vector<std::unique_ptr<Particle>>