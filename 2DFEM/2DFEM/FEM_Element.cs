using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{
    abstract class FEM_Element
    {
        public Node[] Nodes   // 要素内のノード
        {
            get;
            protected set;
        }
        public abstract DenseVector[] NodeStressVector
        {
            get;
            protected set;
        }

        public abstract DenseMatrix makeKeMatrix();

        public abstract void makeStrainVector(DenseVector dispvector);

        public abstract void makeStressVector();

        public abstract void makeNodeStressVector();

        public FEM_Element ShallowCopy()
        {
            return (FEM_Element)MemberwiseClone();
        }
    }
}
