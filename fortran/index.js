const fs = require('fs');

let print_time_count;
let time = 0.0;

let lo_mi, H_mi = 0.0;
let lo_aer, cur_lc_aer = 0.0;

let A1_i = new Array(10).fill(0),
	A2_i = new Array(10).fill(0),
	A3_i = new Array(10).fill(0),
	A4_i = new Array(10).fill(0),
	A1prev_i,
	A2prev_i = new Array(10).fill(0),
	A3prev_i = new Array(10).fill(0),
	A4prev_i = new Array(10).fill(0),
	dA1_i = new Array(10).fill(0),
	dA2_i = new Array(10).fill(0),
	dA3_i = new Array(10).fill(0),
	dA4_i = new Array(10).fill(0),
	A1in_i,
	filter_i = new Array(10).fill(0),
	filter_aer = new Array(10).fill(0),
	Aaccum_i = new Array(10).fill(0);

let lo_i = new Array(10).fill(0),
	ld_i = new Array(10).fill(0),
	lc_i = new Array(10).fill(0),
	cur_lc_i = new Array(10).fill(0),
	H_i = new Array(10).fill(0);

let A1_aer = new Array(7).fill(0),
	A2_aer = new Array(7).fill(0),
	A3_aer = new Array(7).fill(0),
	A4_aer = new Array(7).fill(0),
	A1prev_aer,
	A2prev_aer = new Array(7).fill(0),
	A3prev_aer = new Array(7).fill(0),
	A4prev_aer = new Array(7).fill(0),
	dA1_aer = new Array(7).fill(0),
	dA2_aer = new Array(7).fill(0),
	dA3_aer = new Array(7).fill(0),
	dA4_aer = new Array(7).fill(0),
	A1in_aer,
	Aaccum_aer = new Array(7).fill(0);

let A1in_gas;

let Aout_i  = new Array(10).fill(0).map(
	() => new Array(10).fill(0)),
	A1out_i  = new Array(10).fill(0).map(
		() => new Array(10).fill(0));
let Aout_aer  = new Array(7).fill(0).map(
	() => new Array(10).fill(0));
let Aout_gas  = new Array(6).fill(0).map(
	() => new Array(10).fill(0));

const fileReader = (filename) => {
	const input = {};
	const lines = fs.readFileSync(filename, 'utf-8')
		.split('\n')
		.filter(str => Boolean(str) && str[0] !== '!' && str !== '\\r')
		.map(str => {
			str = str.replace('= ', '');
			str = str.replace('\r', '');
			return str.split( ' ').filter(Boolean);
		});

	lines.forEach(([key, ...values]) => {
		if (key === '/') return;

		if (values.length === 1) {
			input[key] = Number(values[0]);
		} else {
			input[key] = [];
			values.forEach(value => input[key].push(Number(value)));
		}
	});

	return input;
};

let { ld_mi, lc_mi, lc_oi, lc_aer, ly, V1, V2, S, ph, t_s, H_oi, km_mi, km_aer, dt, fin_time, print_time } = fileReader('./fortran/RLS_INP.dat');
let { I1, I2, I3, A1, A2, A3, G1, G2, G3, lr_i, lr_cs, lr_gas } = fileReader('./fortran/RLS_ACT.dat');

const sumVector = (Vec1, Vec2) => Vec1.map((a, i) => a + Vec2[i]);
const sumScalar = (Vec1, scalar) => Vec1.map( a => a + scalar);
const subVector = (Vec1, Vec2) => Vec1.map((a, i) => a - Vec2[i]);
const mulVector = (Vec1, Vec2) => Vec1.map((a, i) => a * Vec2[i]);
const mulScalar = (Vec1, scalar) => Vec1.map( a => a * scalar);
const divVector = (Vec1, Vec2) => Vec1.map( (a, i) => a / Vec2[i]);
const divScalarInv = (scalar, Vec1) => Vec1.map( a => scalar / a);

A1in_i = sumVector(sumVector(I1, I2), I3);
A1in_aer = sumVector(sumVector(A1, A2), A3);
A1in_gas = sumVector(sumVector(G1, G2), G3);

print_time_count = print_time.length;

A1prev_i = [...A1in_i];
A1prev_aer = [...A1in_aer];

ly = ly/(24*3600*100);

lo_mi = km_mi * S / V1;
lo_aer = km_aer * S / V1;

let lr_aer = [...lr_i, lr_cs[0], lr_cs[1]];

if (ph < 6.5) {
	H_mi = 50;
}
else if (ph > 8.5) {
	H_mi = 5000;
}
else {
	H_mi = 50 + (4950 / 2) * (ph-6.5);
}

for (let i = 0; i <= 9; ++i) {
	if (i < 5) {
		lo_i[i] = lo_mi;
		ld_i[i] = ld_mi;
		lc_i[i] = lc_mi;
		H_i[i] = H_mi;
	} else {
		lo_i[i] = 0;
		ld_i[i] = 0;
		lc_i[i] = lc_oi;
		lr_i[i] = lr_i[i - 4];
		H_i[i] = H_oi;
	}
}

while (time < fin_time) {
	time = time + dt;

	if ((time >= t_s) && (time <= t_s + dt)) {
		cur_lc_i = lc_i;
		cur_lc_aer = lc_aer;
	}

	dA1_i = sumVector(
		mulScalar(mulVector(sumScalar(sumVector(sumVector(lo_i, lr_i), cur_lc_i), ly), A1prev_i), -1),
		sumVector(
			mulVector(mulVector(cur_lc_i, divScalarInv(V1, mulScalar(H_i, V2))), A3prev_i),
			mulVector(ld_i, A3prev_i)
		)
	);
	dA2_i = mulScalar(A1prev_i, ly);
	dA3_i = subVector(mulVector(lo_i, A1prev_i), mulVector(sumVector(ld_i, lr_i), A3prev_i));
	dA4_i = subVector(mulVector(cur_lc_i, A1prev_i),mulVector(sumVector(mulVector(cur_lc_i, divScalarInv(V1, mulScalar(H_i, V2))), lr_i), A4prev_i));

	dA1_aer = mulScalar(mulVector(sumScalar(sumScalar(sumScalar(lr_aer, lo_aer), cur_lc_aer), ly), A1prev_aer), -1);
	dA2_aer = mulScalar(A1prev_aer, ly);
	dA3_aer = subVector(mulScalar(A1prev_aer, lo_aer), mulVector(lr_aer, A3prev_aer));
	dA4_aer = subVector(mulScalar(A1prev_aer, cur_lc_aer), mulVector(lr_aer, A4prev_aer));

	A1_i = sumVector(A1prev_i, mulScalar(dA1_i, dt));
	A2_i = sumVector(A2prev_i, mulScalar(dA2_i, dt));
	A3_i = sumVector(A3prev_i, mulScalar(dA3_i, dt));
	A4_i = sumVector(A4prev_i, mulScalar(dA4_i, dt));
	A1_aer = sumVector(A1prev_aer, mulScalar(dA1_aer, dt));
	A2_aer = sumVector(A2prev_aer, mulScalar(dA2_aer, dt));
	A3_aer = sumVector(A3prev_aer, mulScalar(dA3_aer, dt));
	A4_aer = sumVector(A4prev_aer, mulScalar(dA4_aer, dt));

	Aaccum_i = sumVector(Aaccum_i, mulScalar(mulVector(A2_i, sumScalar(mulScalar(filter_i, -1), 1)), dt));
	Aaccum_aer = sumVector(Aaccum_aer, mulScalar(mulVector(A2_aer, sumScalar(mulScalar(filter_aer, -1), 1)), dt));

	for (let i = 0; i < print_time_count; ++i) {
		if (time >= print_time[i]) {
			if ((time - dt) < print_time[i]) {
				Aout_i[i] = mulVector(A2_i, sumScalar(mulScalar(filter_i, -1), 1));
				A1out_i[i] = [...A1_i];
				Aout_aer[i] = mulVector(A2_aer, sumScalar(mulScalar(filter_aer, -1), 1));
				Aout_gas[i] = mulVector(divVector(mulScalar(A1in_gas, ly), sumScalar(lr_gas, ly)),
					sumScalar(mulScalar(mulScalar(sumScalar(lr_gas, ly), (-1 * time)).map(Math.exp), -1), 1));
			}
		}
	}

	A1prev_i = [...A1_i];
	A2prev_i = [...A2_i];
	A3prev_i = [...A3_i];
	A4prev_i = [...A4_i];
	A1prev_aer = [...A1_aer];
	A2prev_aer = [...A2_aer];
	A3prev_aer = [...A3_aer];
	A4prev_aer = [...A4_aer];
}

console.log('TIME', '3600.0', '28800.0', '43200.0', '86400.0', '259200.0');
for (let i = 0; i < 5; ++i) {
	console.log(Aout_i[0][i], Aout_i[1][i], Aout_i[2][i], Aout_i[3][i], Aout_i[4][i]);
}

console.log('\n');
for (let i = 5; i < 10; ++i) {
	console.log(Aout_i[0][i], Aout_i[1][i], Aout_i[2][i], Aout_i[3][i], Aout_i[4][i]);
}

console.log('\n');
for (let i = 0; i < 5; ++i) {
	console.log(Aout_aer[0][i], Aout_aer[1][i], Aout_aer[2][i], Aout_aer[3][i], Aout_aer[4][i]);
}

console.log('\n');
for (let i = 0; i < 5; ++i) {
	console.log(Aout_gas[0][i], Aout_gas[1][i], Aout_gas[2][i], Aout_gas[3][i], Aout_gas[4][i]);
}

// вывод всех данных
/*Aout_i.forEach(row => {
	row.forEach( val => console.log(val));
	console.log('\n');
});
Aout_aer.forEach(row => {
	row.forEach( val => console.log(val));
	console.log('\n');
});
Aout_gas.forEach(row => {
	row.forEach( val => console.log(val));
	console.log('\n');
});*/
