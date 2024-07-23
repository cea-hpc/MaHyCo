// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef QNEWT_STDPERFECTGAS_UTILS_H
#define QNEWT_STDPERFECTGAS_UTILS_H

using namespace Arcane;

namespace Qnewt_stdperfectgas {

template<typename T>
void
resizeObjectFromMask(Arcane::ConstArrayView<bool> conv_mask, T* resizeable)
{
  Integer i_item=0;
  resizeable->erase(
      std::remove_if(
	resizeable->begin(), resizeable->end(),
	[&conv_mask, &i_item]([[maybe_unused]] const typename T::value_type& val)
	{ 
	return conv_mask[i_item++]; 
	}
	),
      resizeable->end()
      );
}

} // Qnewt_stdperfectgas

#endif

