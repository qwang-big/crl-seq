var width = 900, height = 900, margin = 50, rank, highlights=[]
var tooltip = d3.select("#tooltip")
  var chart = d3.select("#chart")
    .append("svg:svg")
    .attr("class", "chart")
    .attr("width", width)
    .attr("height", height);
  var y = d3.scale.linear()
    .domain([1,17888])
    .range([height-margin, margin]);
  var x = d3.scale.linear()
    .domain([1,17888])
    .range([margin,width-margin]);

  chart.selectAll("line.x")
   .data(x.ticks(10))
   .enter().append("svg:line")
   .attr("class", "x")
   .attr("x1", x)
   .attr("x2", x)
   .attr("y1", margin)
   .attr("y2", height - margin)
   .attr("stroke", "#ccc");

  chart.selectAll("line.y")
   .data(y.ticks(10))
   .enter().append("svg:line")
   .attr("class", "y")
   .attr("x1", margin)
   .attr("x2", width - margin)
   .attr("y1", y)
   .attr("y2", y)
   .attr("stroke", "#ccc");

  chart.selectAll("text.xrule")
   .data(x.ticks(10))
   .enter().append("svg:text")
   .attr("class", "xrule")
   .attr("x", x)
   .attr("y", height - margin)
   .attr("dy", 20)
   .attr("text-anchor", "middle")
   .text(String);

 chart.selectAll("text.yrule")
  .data(y.ticks(10))
  .enter().append("svg:text")
  .attr("class", "yrule")
  .attr("x", width - margin+1)
  .attr("y", y)
  .attr("dy", 0)
  .attr("dx", 20)     
  .attr("text-anchor", "middle")
  .text(String);

 /*chart.selectAll("text.yrule2")
  .data(y.ticks(10))
  .enter().append("svg:text")
  .attr("class", "yrule")
  .attr("x", 1)
  .attr("y", y)
  .attr("dy", 0)
  .attr("dx", 20)     
  .attr("text-anchor", "middle")
  .text(String);*/

  chart.append("svg:rect")
  .attr("x", margin)
  .attr("y", margin)      
    .attr("width", width - margin*2)
    .attr("height", height - margin*2)
    .attr("fill-opacity",0)
    .attr('stroke', '#000')
    .attr('stroke-width',1)

    chart.append('text')
    .attr('x',4)
    .attr('y',height/4)
    .text("Low significance <- - - - - - - - - - Promoter only - - - - - - - - - -> High significance")
  .attr("transform", "translate(-210,750) rotate(-90 0 0)")
    .style("font-size", "16px")

    chart.append('text')
    .attr('x',width/4 - margin)
    .attr('y',height-4)
    .text("Low significance <- - - - - - - - - - Promoter + Enhancer - - - - - - - - - -> High significance")
    .style("font-size", "16px")

function highlight(i){
	chart.selectAll(".bar").remove()
	d3.csv(csv[i], function(error,data){
		highlights=[]
		for(i in data) highlights.push(data[i].name)
		plot(rank)
	})
}

function plot(data){
  var bar = chart.selectAll(".bar")
    .data(data)
    .enter().append("g")
    .attr("class", "bar");

  bar.append("svg:rect")
    .on('click', function(d,i) {alert(d.promenh)})
    .attr("x", function(d) { return x(d.promenh) })
    .attr("y", function(d) {return y(d.promonly)})      
    .attr("height", function(d) { return highlights.indexOf(d.name)<0 ? 3 : 7})
    .attr("width", function(d) { return highlights.indexOf(d.name)<0 ? 3 : 7})
    .attr("fill",function(d) { return highlights.indexOf(d.name)<0 ? "#888" : "#f00"})
  .on("mouseover", function(d){
	tooltip.html(d.name)
	return tooltip.style("visibility", "visible");
      })
      .on("mousemove", function(){
           return tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px")
      })
      .on("mouseout", function(){
           return tooltip.style("visibility", "hidden")
      })
}

d3.csv(csv[0],function(error,data){
	prom =[], rank =[]
	for(i in data) {
		prom[data[i].promonly]=data.length-i
	}
	for(i in data) {
		rank.push({"name": data[i].promenh, "promenh": data.length-i, "promonly": prom[data[i].promenh]})
	}
	plot(rank);
})
