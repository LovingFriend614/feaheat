var N = document.getElementById("Nin").value;
var canvas = document.getElementById("elementSquare");
var ctx = canvas.getContext("2d");

function linesDraw(Nin) {
	ctx.clearRect(0,0,500,500);
	var N = document.getElementById("Nin").value;
	var changeVal = 500 / N;
	for (i = 0; i < N + 1; i++) {
		var addVal = changeVal * i;
		ctx.beginPath();
		ctx.moveTo(0,addVal);
		ctx.lineTo(500,-500 + addVal);
		ctx.moveTo(addVal,0);
		ctx.lineTo(addVal,500);
		ctx.moveTo(0,addVal);
		ctx.lineTo(500,addVal);
		ctx.closePath();
		ctx.lineWidth = 3;
		ctx.stroke();
	}
}
