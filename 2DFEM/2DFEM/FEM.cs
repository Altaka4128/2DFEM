using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{
    abstract class FEM
    {
        public abstract List<FEM_Element> Elems   // 要素のリスト
        {
            get;
            protected set;
        }
        public DenseVector DispVector   // 変位の境界条件
        {
            get;
            protected set;
        }
        public DenseVector ForceVector   // 荷重の境界条件
        {
            get;
            protected set;
        }
        public List<bool> Constraint   // 拘束の境界条件
        {
            get;
            protected set;
        }

        public abstract void Analysis();
    }
}
