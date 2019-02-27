NAME=midas
TAG=$$(date +%Y%m%d%H%M%S)
REGISTRY=sunitjain/$(NAME)

all: build push run

build:
	docker build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile .
	docker tag $(NAME):latest $(REGISTRY):latest
	docker tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker push $(REGISTRY):$(TAG)
	docker push $(REGISTRY):latest

run:
	docker run $(NAME):latest
