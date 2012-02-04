template class Large_Matrix
{
private:
	unsigned nr_rows;
	unsigned nr_columns;
	T * data;
public:
	Large_Matrix(unsigned nr_rows_,unsigned nr_columns_) :
		nr_rows(nr_rows_),
		nr_columns(nr_columns_),
		data(new T[nr_rows_*nr_columns_])
	{}
	~Large_Matrix()
	{
		delete [] data;
	}
	T * operator[](unsigned row)
	{
		return &data[row*nr_columns];
	}
	const T * operator[](unsigned row) const
	{
		return &data[row*nr_columns];
	}
	static bool multiply(Large_Matrix * answer,
						 const Large_Matrix & lhs,
						 const Large_Matrix & rhs)
	{
		if (lhs.nr_columns != rhs.nr_rows)
			return false;
		if (answer == &lhs)
		{
			Large_Matrix temp(lhs);
			return multiply(answer,temp,rhs);
		}
		if (answer == &rhs)
		{
			Large_Matrix temp(rhs);
			return multiply(answer,lhs,temp);
		}
		unsigned nr_elements = lhs.nr_rows * rhs.nr_columns;
		if (answer->nr_rows * answer->nr_columns != nr_elements)
		{
			delete [] answer->data;
			answer->data = new T[lhs.nr_rows*rhs.nr_columns];
		}
		answer->nr_rows = lhs.nr_rows;
		answer->nr_columns = rhs.nr_columns;

		{
			for (unsigned i = 0; i data[i] = 0;
				 }

			unsigned i0,j0,k0,i,j,k,imax,jmax,kmax,jmax1;
			const unsigned IBLOCK=100;
			const unsigned JBLOCK=128;
			const unsigned KBLOCK=40;
			for (i0 = 0; i0 lhs.nr_rows)
				imax = lhs.nr_rows;
			for (k0 = 0; k0 lhs.nr_columns)
				kmax = lhs.nr_columns;
			for (j0 = 0; j0 rhs.nr_columns)
				jmax = rhs.nr_columns;
			jmax1 = jmax & ~15;
			for (i = i0; i <imax;i++,lrow += lhs.nr_columns,arow += rhs.nr_columns)
			{
				const T * rrow = rhs[k0];
				for (k = k0; k < kmax;k++,rrow += rhs.nr_columns)
				{
					T v = lrow[k];
					for (j = j0; j < jmax1;j+=16)
					{
						arow[j] += v*rrow[j];
						arow[j+1] += v*rrow[j+1];
						arow[j+2] += v*rrow[j+2];
						arow[j+3] += v*rrow[j+3];
						arow[j+4] += v*rrow[j+4];
						arow[j+5] += v*rrow[j+5];
						arow[j+6] += v*rrow[j+6];
						arow[j+7] += v*rrow[j+7];
						arow[j+8] += v*rrow[j+8];
						arow[j+9] += v*rrow[j+9];
						arow[j+10] += v*rrow[j+10];
						arow[j+11] += v*rrow[j+11];
						arow[j+12] += v*rrow[j+12];
						arow[j+13] += v*rrow[j+13];
						arow[j+14] += v*rrow[j+14];
						arow[j+15] += v*rrow[j+15];
					}
					for (j = jmax1; j < jmax;j++)
						arow[j] += v*rrow[j];
				}
			}
		}
	}
}
	return true;
}
};
