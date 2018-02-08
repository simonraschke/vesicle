/*  
*   Copyright 2017-2018 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#pragma once

#include "translator.hpp"
#include "particles/particle_mobile.hpp"
#include "vesicleIO/gro_reader.hpp"
#include <type_traits>
#include <string>




class AnisotropicCoordsTranslatorGro
    : public Translator<Particle::cartesian>
{
public:
    virtual Particle::cartesian operator()(std::string,std::string) const override;
};



class AnisotropicOrientationTranslatorGro
    : public Translator<Particle::cartesian>
{
public:
    virtual Particle::cartesian operator()(std::string,std::string) const override;
};



class AnisotropicParticleTranslatorGro
    : public Translator<Particle*>
{
public:
    virtual Particle* operator()(std::string,std::string) const override;
};



class IsotropicParticleTranslatorGro
    : public Translator<Particle*>
{
public:
    virtual Particle* operator()(std::string) const override;
};