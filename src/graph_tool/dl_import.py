#! /usr/bin/env python
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007 Tiago de Paula Peixoto <tiago@forked.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
try:
    from dl import RTLD_LAZY, RTLD_NOW, RTLD_GLOBAL
except ImportError:
    RTLD_LAZY = 1
    RTLD_NOW = 2
    RTLD_GLOBAL = 256

all = ["dl_import"]

def dl_import(import_expr):
    """Import module according to import_expr, but with RTLD_GLOBAL enabled."""
    # we need to get the locals and globals of the _calling_ function. Thus, we
    # need to go deeper into the call stack
    call_frame = sys._getframe(1)
    local_dict = call_frame.f_locals
    global_dict = call_frame.f_globals

    # RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and friends to
    # work properly across DSO boundaries. See http://gcc.gnu.org/faq.html#dso

    # The "except" is because the dl module raises a system error on ia64 and x86_64
    # systems because "int" and addresses are different sizes.
    orig_dlopen_flags = sys.getdlopenflags()

    sys.setdlopenflags(RTLD_LAZY|RTLD_GLOBAL)

    exec import_expr in local_dict, global_dict

    sys.setdlopenflags(orig_dlopen_flags) # reset it to normal case to avoid
                                          # unnecessary symbol collision
