.PHONY: build

build:
	docker build -t mixaill76/poly-vectorization:latest .

run:
	docker run -it --rm \
		-e DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
		-v $XAUTHORITY:/root/.Xauthority:ro \
		-v ./sample_inputs:/app/sample_inputs \
		mixaill76/poly-vectorization:latest