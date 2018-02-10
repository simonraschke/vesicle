#include "cluster.hpp"



void Cluster::add(Particle& particle)
{
    mutex_type::scoped_lock lock(mutex);
    members.emplace_back(std::addressof(particle));
}



bool Cluster::contains(const Particle* other) const
{
    // auto particle_oberver = enhance::make_observer<Particle>(particle);
    return std::none_of(begin(),end(),[&](const MemberType& member){ return member == other; });
}



bool Cluster::contains(const ParticleSimple& other) const
{
    // auto particle_oberver = enhance::make_observer<Particle>(particle);
    return std::none_of(begin(),end(),[&](const MemberType& member){ return member == other; });
}



void Cluster::add_if_new(Particle& particle)
{
    if( !contains(std::addressof(particle)) ) 
    {
        add(particle);
    }
}



// void Cluster::expandNeighbourhood(MemberNeighbourPair& pair, Particle* particle)
// {
//     mutex_type::scoped_lock lock(mutex);
//     pair.second.emplace_back(enhance::make_observer(particle));
// }



void Cluster::remove(Particle* other)
{
    mutex_type::scoped_lock lock(mutex);
    // auto particle_oberver = enhance::make_observer<Particle>(particle);
    members.erase(std::remove(begin(), end(), other), end());
}



void Cluster::remove(const ParticleSimple& other)
{
    mutex_type::scoped_lock lock(mutex);
    // auto particle_oberver = enhance::make_observer<Particle>(particle);
    members.erase(std::remove(begin(), end(), other), end());
}



Cluster::MemberList::iterator Cluster::begin()
{
    return std::begin(members);
}



Cluster::MemberList::iterator Cluster::end()
{
    return std::end(members);
}



Cluster::MemberList::const_iterator Cluster::begin() const
{
    return std::cbegin(members);
}



Cluster::MemberList::const_iterator Cluster::end() const
{
    return std::cend(members);
}