/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "szv_grid.h"
#include "brick.h"

szv_grid::szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr)
    : m_data(gd[0].n, gd[1].n, gd[2].n),
      m_init(gd[0].begin, gd[1].begin, gd[2].begin),
      m_end(gd[0].end, gd[1].end, gd[2].end),
      m_range(m_end - m_init)
{
	const sz nat = num_atom_types(m.atom_typing_used());

	szv relevant_indexes;
	VINA_FOR_IN(i, m.grid_atoms) {
		const atom& a = m.grid_atoms[i];
		if(a.get(m.atom_typing_used()) < nat && brick_distance_sqr(m_init, m_end, a.coords) < cutoff_sqr)
			relevant_indexes.push_back(i);
	}

	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		VINA_FOR_IN(ri, relevant_indexes) {
			const sz i = relevant_indexes[ri];
			const atom& a = m.grid_atoms[i];
			if(brick_distance_sqr(index_to_coord(x, y, z), index_to_coord(x+1, y+1, z+1), a.coords) < cutoff_sqr)
				m_data(x, y, z).push_back(i);
		}
	}
}

fl szv_grid::average_num_possibilities() const {
	sz counter = 0;
	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		counter += m_data(x, y, z).size();
	}
	return fl(counter) / (m_data.dim0() * m_data.dim1() * m_data.dim2());
}

const szv& szv_grid::possibilities(const vec& coords) const {
	boost::array<sz, 3> index;
    VINA_FOR(i, index.static_size) {
		assert(coords[i] + epsilon_fl >= m_init[i]);
		assert(coords[i] <= m_init[i] + m_range[i] + epsilon_fl);
		const fl tmp = (coords[i] - m_init[i]) * m_data.dim(i) / m_range[i];
		index[i] = fl_to_sz(tmp, m_data.dim(i) - 1);
	}
	return m_data(index[0], index[1], index[2]);
}

vec szv_grid::index_to_coord(sz i, sz j, sz k) const {
    assert(m_init.size() == m_range.size() && m_range.size() == 3);
	return vec(m_init[0] + m_range[0] * i / m_data.dim0(),
               m_init[1] + m_range[1] * j / m_data.dim1(),
               m_init[2] + m_range[2] * k / m_data.dim2());
}

grid_dims szv_grid_dims(const grid_dims& gd) {
	grid_dims tmp;
    VINA_FOR(i, tmp.static_size) {
		tmp[i].begin = gd[i].begin;
		tmp[i].end   = gd[i].end;
		const fl n_fl = (gd[i].end - gd[i].begin) / 3; // 3A preferred size
		const int n_int = int(n_fl);
		tmp[i].n = (n_int < 1) ?  1 : sz(n_int);
	}
	return tmp;
}
