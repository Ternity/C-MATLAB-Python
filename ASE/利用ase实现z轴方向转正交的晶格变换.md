参考pymatgen的做法，源代码为：

```
def get_orthogonal_c_slab(self) -> Self:
    """Generate a Slab where the normal (c lattice vector) is forced to be orthogonal to the surface a and b lattice vectors.
    **Note that this breaks inherent symmetries in the slab.**
    It should be pointed out that orthogonality is not required to get good surface energies, but it can be useful in cases where the slabs are subsequently used for postprocessing of some kind, e.g. generating grain boundaries or interfaces. """
    a, b, c = self.lattice.matrix
    _new_c = np.cross(a, b)             # 求(a,b)叉积，得垂直于(a,b)的向量
    _new_c /= np.linalg.norm(_new_c)    # 先求模(范数), 再做归一化,得单位向量
    new_c = np.dot(c, _new_c) * _new_c # 先求向量c在单位向量投影长度, ↓
    new_latt = Lattice([a, b, new_c]) # 再赋值到单位向量方向上
 
    return type(self)(
        lattice=new_latt,
        species=self.species_and_occu,
        coords=self.cart_coords,
        miller_index=self.miller_index,
        oriented_unit_cell=self.oriented_unit_cell,
        shift=self.shift,
        scale_factor=self.scale_factor,
        coords_are_cartesian=True,
        energy=self.energy,
        reorient_lattice=self.reorient_lattice,
        site_properties=self.site_properties,
    )
```

pymatgen在实现这个函数时候，返回的类里面使用了coords\_are\_cartesian\=True方法,这个方法(默认False)会将原子位置和新的盒子矩阵的逆矩阵(np.linalg.inv)做一次点乘dot,从而更新原子位置。



那么在ase中，虽然没有直接实现使α和β\=90°的method，但是可以参考pymatgen实现：

```
import numpy as np

old_crystal = old_Atoms_class
a, b, c = old_crystal.cell
_new_c = np.cross(a, b)            # 求(a,b)叉积，得垂直于(a,b)的向量
_new_c /= np.linalg.norm(_new_c)   # 先求模(范数), 再做归一化,得单位向量
new_c = np.dot(c, _new_c) * _new_c # 先求向量c在单位向量投影长度, 再×单位向量方向

# method 1      
ols_crystal.set_cell([a, b, new_c], scale_atoms=True) # scale_atoms根据盒子变坐标

# method 2
old_crystal.set_cell([a, b, new_c])
M = np.linalg.solve(old_Atoms_class.cell.complete(), old_crystal.cell.complete())   # 求A到B的转换矩阵M
old_crystal.positions[:] = np.dot(old_crystal.positions, M) # [:] 切片操作,就地修改, 直接赋值
```

