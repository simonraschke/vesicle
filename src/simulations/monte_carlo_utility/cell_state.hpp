/*  
*   Copyright 2017 Simon Raschke
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

#include <atomic>


class CellState
{
public:
    enum STATE { UNDEFINED, IDLE, BLOCKED, FINISHED };

    // CellState& operator=(const CellState&);
    CellState& operator=(const STATE&);
    bool operator==(const CellState&) const;
    bool operator==(STATE) const;
    bool operator!=(const CellState&) const;
    bool operator!=(STATE) const;
private:
    std::atomic<STATE> state {UNDEFINED};
};



// inline CellState& CellState::operator=(const CellState& s)
// {
//     state.store(s.state);
//     return *this;
// }



inline CellState& CellState::operator=(const STATE& s)
{
    state.store(s);
    return *this;
}



inline bool CellState::operator==(const CellState& other) const
{
    return state.load() == other.state.load();
}



inline bool CellState::operator==(STATE other) const
{
    return state.load() == other;
}



inline bool CellState::operator!=(const CellState& other) const
{
    return state.load() != other.state.load();
}



inline bool CellState::operator!=(STATE other) const
{
    return state.load() != other;
}
