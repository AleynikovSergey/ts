const fs = require('fs');
const path = require('path');

const SECONDS_IN_HOUR = 3600,
	  SECONDS_IN_MINUTE  = 60;

const fileFilter = filename => fs.readFileSync(filename, 'utf-8').split('\n')
	.filter(str => Boolean(str) && str[0] !== '#' && str !== '\\r')
	.map(str => {
		str = str.replace('= ', '');
		str = str.replace('\r', '');
		return str.split( ' ').filter(Boolean);
	});

const fileReader = filename => {
	const input = {};
	const lines = fileFilter(filename);

	lines.forEach(([key, value]) => {
		if (Number.isNaN(Number(value))) {
			input[key] = value;
		} else {
			input[key] = Number(value);
		}
	});

	return input;
};

const stepsReader = dirname => {
	const steps = [];

	fs.readdirSync(dirname).forEach(
		filename => steps.push(fileReader(path.join(dirname, filename))),
	);

	return steps;
};

const timeToSecond = time => {
	const [hours, minutes] = time.split(':');

	return Number(hours) * SECONDS_IN_HOUR + Number(minutes) * SECONDS_IN_MINUTE;
};

const logger = (key, value) => console.log(`${key}: ${value.toExponential(2)}`);


const CGas = fileReader('./diploma/Gas_concentration.txt'),
	  CWater = fileReader('./diploma/Water_concentration.txt');

const steps = stepsReader('./diploma/steps');

let mATimeBegin = 0
let mBTimeBegin = 0;
let mCTTimeBegin = 0;

steps.forEach((step, index) => {
	const { TimeBegin, TimeEnd, VThai, VWaterEnd, SThai, COWall } = step;
	const deltaT = timeToSecond(TimeEnd) - timeToSecond(TimeBegin);
	let VWaterBegin = VWaterEnd;

	if (step.VWaterBegin) {
		VWaterBegin = step.VWaterBegin;
	}

	if (index < 2) {
		mATimeBegin = VThai * CGas[TimeBegin];
		mBTimeBegin = VWaterEnd * CWater[TimeBegin];
		mCTTimeBegin = 0;

		if (COWall) {
			mCTTimeBegin = COWall * SThai;
		}
	}

	let K0 = -1 * Math.log(CGas[TimeEnd] / CGas[TimeBegin]) / deltaT;
	let mATimeEnd = mATimeBegin * Math.exp( -1 * K0 * deltaT);
	let K1 = (CWater[TimeEnd] * VWaterEnd - CWater[TimeBegin] * VWaterBegin) / (deltaT * mATimeBegin);
	let mBTimeEnd = mBTimeBegin + (K1 * mATimeBegin * (1 - Math.exp(-1 * K0 * deltaT))) / K0;

	let K2 = K0 - K1;

	let mCTCalcTimeEnd = mCTTimeBegin + (K2 * mATimeBegin) / K0 * (1 - Math.exp(-1 * K0 * deltaT));

	console.log(`Проверяем сходимость фазы ${index + 1}`);

	logger('K0', K0);
	logger('K1', K1);
	logger('K2', K2);
	logger('Масса йода в газе:', mATimeEnd);
	logger('Масса йода в воде:', mBTimeEnd);
	logger('Масса йода на стенке:', mCTCalcTimeEnd);
	console.log(Math.round(mATimeEnd + mBTimeEnd + mCTCalcTimeEnd));
	console.log(mATimeEnd + mBTimeEnd + mCTCalcTimeEnd);

	if (index > 0) {
		mATimeBegin = mATimeEnd;
		mBTimeBegin = mBTimeEnd;
		mCTTimeBegin = mCTCalcTimeEnd;
	}
});
